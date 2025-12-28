#!/usr/bin/env julia
#
# @file mfem_cuda_block.jl
# @brief CUDA implementation of 2D MFEM with one block per coarse element
#
# Each CUDA block processes one coarse element:
# - Threads cooperatively assemble fine elements
# - Block-cooperative PCG solves for interior basis functions
# - Block reductions compute coarse element matrices
#
# Uses hybrid shared/global memory approach:
# - Shared memory: PCG workspace (r, z, p, Ap, diag) + reduction scratch
# - Global memory: Sparse matrix values, phi, rhs, workspace arrays
# - Pre-computed: Sparsity patterns (same for all elements with same ratio)
#

using CUDA
using LinearAlgebra
using SparseArrays
using Printf

# Load shared FEM base kernels (for mesh connectivity and coloring)
include("FEMBase.jl")
using .FEMBase

# ============================================================================
# Constants and Configuration
# ============================================================================

# PCG parameters
const PCG_TOL = 1e-12
const PCG_MAXITER = 1000

# Gauss Quadrature (2x2 rule)
const GAUSS_PT = 1.0 / sqrt(3.0)
const GAUSS_PTS_XI = (-GAUSS_PT, GAUSS_PT, GAUSS_PT, -GAUSS_PT)
const GAUSS_PTS_ETA = (-GAUSS_PT, -GAUSS_PT, GAUSS_PT, GAUSS_PT)
const GAUSS_WTS = (1.0, 1.0, 1.0, 1.0)

# ============================================================================
# Device-side Q1 Shape Functions
# ============================================================================

@inline function shape_functions_device(xi::Float64, eta::Float64)
    N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
    N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
    N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
    N4 = 0.25 * (1.0 - xi) * (1.0 + eta)
    return (N1, N2, N3, N4)
end

@inline function shape_gradients_device(xi::Float64, eta::Float64)
    dN1_dxi = -0.25 * (1.0 - eta)
    dN2_dxi =  0.25 * (1.0 - eta)
    dN3_dxi =  0.25 * (1.0 + eta)
    dN4_dxi = -0.25 * (1.0 + eta)
    dN1_deta = -0.25 * (1.0 - xi)
    dN2_deta = -0.25 * (1.0 + xi)
    dN3_deta =  0.25 * (1.0 + xi)
    dN4_deta =  0.25 * (1.0 - xi)
    return (dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN1_deta, dN2_deta, dN3_deta, dN4_deta)
end

# ============================================================================
# Block-Level Reduction Primitives
# ============================================================================

@inline function block_reduce_sum(g, shmem_reduce, val::Float64, blocksize::Int32)
    tid = threadIdx().x
    @inbounds shmem_reduce[tid] = val
    CG.sync(g)

    s = blocksize ÷ Int32(2)
    while s > Int32(0)
        if tid <= s
            @inbounds shmem_reduce[tid] += shmem_reduce[tid + s]
        end
        CG.sync(g)
        s ÷= Int32(2)
    end

    @inbounds result = shmem_reduce[1]
    CG.sync(g)
    return result
end

@inline function block_dot_product_global(g, shmem_reduce, x, y, n::Int32, blocksize::Int32)
    tid = threadIdx().x
    local_sum = 0.0
    i = tid
    while i <= n
        @inbounds local_sum += x[i] * y[i]
        i += blocksize
    end
    return block_reduce_sum(g, shmem_reduce, local_sum, blocksize)
end

@inline function block_spmv_global!(g, y, n::Int32, colptr, rowidx, nzval, x, blocksize::Int32)
    tid = threadIdx().x
    col = tid
    while col <= n
        @inbounds col_start = colptr[col]
        @inbounds col_end = colptr[col + 1] - Int32(1)
        sum_val = 0.0
        for k = col_start:col_end
            @inbounds row = rowidx[k]
            @inbounds sum_val += nzval[k] * x[row]
        end
        @inbounds y[col] = sum_val
        col += blocksize
    end
    CG.sync(g)
end

# ============================================================================
# Block-Cooperative PCG Solver (Hybrid Memory Version)
# ============================================================================

"""
Block-cooperative PCG with Jacobi preconditioner.
Uses shared memory for workspace (r, z, p, Ap) and global memory for matrix values.
"""
@inline function block_pcg_jacobi_hybrid!(g, x, n::Int32,
                                           colptr, rowidx, nzval,  # Global memory
                                           b,                      # Global memory
                                           diag_ii,               # Shared memory
                                           r, z, p, Ap,           # Shared memory
                                           shmem_reduce,          # Shared memory
                                           blocksize::Int32,
                                           tol::Float64, maxiter::Int32)
    tid = threadIdx().x

    # Compute initial residual: r = b - A*x (r in shared, x in global)
    # First compute A*x into r
    col = tid
    while col <= n
        @inbounds col_start = colptr[col]
        @inbounds col_end = colptr[col + 1] - Int32(1)
        sum_val = 0.0
        for k = col_start:col_end
            @inbounds row = rowidx[k]
            @inbounds sum_val += nzval[k] * x[row]
        end
        @inbounds r[col] = b[col] - sum_val
        col += blocksize
    end
    CG.sync(g)

    # Apply Jacobi preconditioner: z = D^{-1} * r, p = z
    i = tid
    while i <= n
        @inbounds d = diag_ii[i]
        @inbounds z_val = (d != 0.0) ? r[i] / d : r[i]
        @inbounds z[i] = z_val
        @inbounds p[i] = z_val
        i += blocksize
    end
    CG.sync(g)

    # Compute initial rz = r'*z (both in shared memory)
    local_sum = 0.0
    i = tid
    while i <= n
        @inbounds local_sum += r[i] * z[i]
        i += blocksize
    end
    rz_old = block_reduce_sum(g, shmem_reduce, local_sum, blocksize)

    # PCG iteration
    for iter = Int32(1):maxiter
        # Ap = A * p (p in shared, nzval in global, Ap in shared)
        col = tid
        while col <= n
            @inbounds col_start = colptr[col]
            @inbounds col_end = colptr[col + 1] - Int32(1)
            sum_val = 0.0
            for k = col_start:col_end
                @inbounds row = rowidx[k]
                @inbounds sum_val += nzval[k] * p[row]
            end
            @inbounds Ap[col] = sum_val
            col += blocksize
        end
        CG.sync(g)

        # pAp = p' * Ap
        local_sum = 0.0
        i = tid
        while i <= n
            @inbounds local_sum += p[i] * Ap[i]
            i += blocksize
        end
        pAp = block_reduce_sum(g, shmem_reduce, local_sum, blocksize)

        if pAp <= 0.0
            return -iter
        end

        alpha = rz_old / pAp

        # x += alpha * p (x in global, p in shared), r -= alpha * Ap
        i = tid
        while i <= n
            @inbounds x[i] += alpha * p[i]
            @inbounds r[i] -= alpha * Ap[i]
            i += blocksize
        end
        CG.sync(g)

        # Check convergence: ||r||^2
        local_sum = 0.0
        i = tid
        while i <= n
            @inbounds local_sum += r[i] * r[i]
            i += blocksize
        end
        rnorm_sq = block_reduce_sum(g, shmem_reduce, local_sum, blocksize)

        if sqrt(rnorm_sq) < tol
            return iter
        end

        # z = D^{-1} * r
        i = tid
        while i <= n
            @inbounds d = diag_ii[i]
            @inbounds z[i] = (d != 0.0) ? r[i] / d : r[i]
            i += blocksize
        end
        CG.sync(g)

        # rz_new = r' * z
        local_sum = 0.0
        i = tid
        while i <= n
            @inbounds local_sum += r[i] * z[i]
            i += blocksize
        end
        rz_new = block_reduce_sum(g, shmem_reduce, local_sum, blocksize)

        beta = rz_new / rz_old

        # p = z + beta * p
        i = tid
        while i <= n
            @inbounds p[i] = z[i] + beta * p[i]
            i += blocksize
        end
        CG.sync(g)

        rz_old = rz_new
    end

    return -maxiter
end

# ============================================================================
# Fine Element Assembly (Device Functions)
# ============================================================================

@inline function assemble_fine_element_device(x1::Float64, y1::Float64, x2::Float64, y2::Float64,
                                               x3::Float64, y3::Float64, x4::Float64, y4::Float64,
                                               use_varying_coeffs::Bool)
    Ke11 = 0.0; Ke12 = 0.0; Ke13 = 0.0; Ke14 = 0.0
    Ke21 = 0.0; Ke22 = 0.0; Ke23 = 0.0; Ke24 = 0.0
    Ke31 = 0.0; Ke32 = 0.0; Ke33 = 0.0; Ke34 = 0.0
    Ke41 = 0.0; Ke42 = 0.0; Ke43 = 0.0; Ke44 = 0.0
    fe1 = 0.0; fe2 = 0.0; fe3 = 0.0; fe4 = 0.0

    for qp = 1:4
        xi = GAUSS_PTS_XI[qp]
        eta = GAUSS_PTS_ETA[qp]
        w = GAUSS_WTS[qp]

        N1, N2, N3, N4 = shape_functions_device(xi, eta)
        dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN1_deta, dN2_deta, dN3_deta, dN4_deta = shape_gradients_device(xi, eta)

        J11 = dN1_dxi * x1 + dN2_dxi * x2 + dN3_dxi * x3 + dN4_dxi * x4
        J12 = dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3 + dN4_deta * x4
        J21 = dN1_dxi * y1 + dN2_dxi * y2 + dN3_dxi * y3 + dN4_dxi * y4
        J22 = dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3 + dN4_deta * y4

        detJ = J11 * J22 - J12 * J21
        invdetJ = 1.0 / detJ

        invJ11 =  J22 * invdetJ
        invJ12 = -J12 * invdetJ
        invJ21 = -J21 * invdetJ
        invJ22 =  J11 * invdetJ

        dN1_dx = dN1_dxi * invJ11 + dN1_deta * invJ21
        dN1_dy = dN1_dxi * invJ12 + dN1_deta * invJ22
        dN2_dx = dN2_dxi * invJ11 + dN2_deta * invJ21
        dN2_dy = dN2_dxi * invJ12 + dN2_deta * invJ22
        dN3_dx = dN3_dxi * invJ11 + dN3_deta * invJ21
        dN3_dy = dN3_dxi * invJ12 + dN3_deta * invJ22
        dN4_dx = dN4_dxi * invJ11 + dN4_deta * invJ21
        dN4_dy = dN4_dxi * invJ12 + dN4_deta * invJ22

        # Compute physical coordinates at quadrature point
        xq = N1 * x1 + N2 * x2 + N3 * x3 + N4 * x4
        yq = N1 * y1 + N2 * y2 + N3 * y3 + N4 * y4

        # Evaluate coefficients at quadrature point (matching CPU version)
        if use_varying_coeffs
            ax_val = 1.0 + 100.0 * cos(150.0 * xq)^2 * sin(150.0 * yq)^2
            ay_val = ax_val
            f_val = sin(xq) * sin(yq)
        else
            ax_val = 1.0
            ay_val = 1.0
            f_val = 2.0 * 3.141592653589793^2 * sin(3.141592653589793 * xq) * sin(3.141592653589793 * yq)
        end

        detJ_w = detJ * w
        ax_detJ_w = ax_val * detJ_w
        ay_detJ_w = ay_val * detJ_w
        f_detJ_w = f_val * detJ_w

        Ke11 += ax_detJ_w * dN1_dx * dN1_dx + ay_detJ_w * dN1_dy * dN1_dy
        Ke12 += ax_detJ_w * dN1_dx * dN2_dx + ay_detJ_w * dN1_dy * dN2_dy
        Ke13 += ax_detJ_w * dN1_dx * dN3_dx + ay_detJ_w * dN1_dy * dN3_dy
        Ke14 += ax_detJ_w * dN1_dx * dN4_dx + ay_detJ_w * dN1_dy * dN4_dy

        Ke21 += ax_detJ_w * dN2_dx * dN1_dx + ay_detJ_w * dN2_dy * dN1_dy
        Ke22 += ax_detJ_w * dN2_dx * dN2_dx + ay_detJ_w * dN2_dy * dN2_dy
        Ke23 += ax_detJ_w * dN2_dx * dN3_dx + ay_detJ_w * dN2_dy * dN3_dy
        Ke24 += ax_detJ_w * dN2_dx * dN4_dx + ay_detJ_w * dN2_dy * dN4_dy

        Ke31 += ax_detJ_w * dN3_dx * dN1_dx + ay_detJ_w * dN3_dy * dN1_dy
        Ke32 += ax_detJ_w * dN3_dx * dN2_dx + ay_detJ_w * dN3_dy * dN2_dy
        Ke33 += ax_detJ_w * dN3_dx * dN3_dx + ay_detJ_w * dN3_dy * dN3_dy
        Ke34 += ax_detJ_w * dN3_dx * dN4_dx + ay_detJ_w * dN3_dy * dN4_dy

        Ke41 += ax_detJ_w * dN4_dx * dN1_dx + ay_detJ_w * dN4_dy * dN1_dy
        Ke42 += ax_detJ_w * dN4_dx * dN2_dx + ay_detJ_w * dN4_dy * dN2_dy
        Ke43 += ax_detJ_w * dN4_dx * dN3_dx + ay_detJ_w * dN4_dy * dN3_dy
        Ke44 += ax_detJ_w * dN4_dx * dN4_dx + ay_detJ_w * dN4_dy * dN4_dy

        fe1 += N1 * f_detJ_w
        fe2 += N2 * f_detJ_w
        fe3 += N3 * f_detJ_w
        fe4 += N4 * f_detJ_w
    end

    return (Ke11, Ke12, Ke13, Ke14,
            Ke21, Ke22, Ke23, Ke24,
            Ke31, Ke32, Ke33, Ke34,
            Ke41, Ke42, Ke43, Ke44,
            fe1, fe2, fe3, fe4)
end

@inline function find_matrix_position_device(rowptr, colidx, row::Int32, col::Int32)
    @inbounds left = rowptr[row]
    @inbounds right = rowptr[row + 1] - Int32(1)

    while left <= right
        mid = (left + right) >> 1
        @inbounds col_mid = colidx[mid]
        if col_mid == col
            return mid
        elseif col_mid < col
            left = mid + Int32(1)
        else
            right = mid - Int32(1)
        end
    end
    return Int32(-1)
end

# ============================================================================
# Main MFEM Element Kernel - Hybrid Memory Version
# ============================================================================

"""
MFEM element assembly kernel - one block per coarse element (hybrid memory version)

Shared memory layout (only PCG workspace):
- diag_ii: nfree Float64
- r, z, p, Ap: 4 * nfree Float64
- shmem_reduce: blocksize Float64

Global memory (per-element workspaces, indexed by iel):
- nzval_ii: max_nnz_ii per element
- matValues_b: max_nnz_b per element
- phi: numNodes * 4 per element
- rhs_fine: numNodes per element
- btmp: nfree * 3 per element
- utmp: nfree * 3 per element

Pre-computed (shared across all elements):
- Sparsity patterns: colptr_ii, rowidx_ii, matRowPtr_b, matColIdx_b
- DOF mappings: globalToFree, freeToGlobal, globalToBoundary, boundaryToGlobal
- Boundary phi values: phi_boundary
"""
function mfem_element_kernel_hybrid!(
    # Outputs
    Ke_out::CuDeviceMatrix{Float64},
    fe_out::CuDeviceMatrix{Float64},
    phi_out::CuDeviceArray{Float64, 3},
    # Mesh data
    vertex_x::CuDeviceVector{Float64},
    vertex_y::CuDeviceVector{Float64},
    cell_to_node::CuDeviceMatrix{Int32},
    # Pre-computed sparsity (shared across elements)
    colptr_ii::CuDeviceVector{Int32},
    rowidx_ii::CuDeviceVector{Int32},
    matRowPtr_b::CuDeviceVector{Int32},
    matColIdx_b::CuDeviceVector{Int32},
    # Pre-computed DOF mappings (shared across elements)
    globalToFree::CuDeviceVector{Int32},
    freeToGlobal::CuDeviceVector{Int32},
    globalToBoundary::CuDeviceVector{Int32},
    boundaryToGlobal::CuDeviceVector{Int32},
    # Pre-computed boundary phi values
    phi_boundary::CuDeviceMatrix{Float64},  # nboundary × 4
    # Per-element workspaces (global memory)
    nzval_ii_all::CuDeviceMatrix{Float64},      # nel × max_nnz_ii
    matValues_b_all::CuDeviceMatrix{Float64},   # nel × max_nnz_b
    phi_all::CuDeviceArray{Float64, 3},         # nel × numNodes × 4
    rhs_fine_all::CuDeviceMatrix{Float64},      # nel × numNodes
    btmp_all::CuDeviceArray{Float64, 3},        # nel × nfree × 3
    utmp_all::CuDeviceArray{Float64, 3},        # nel × nfree × 3
    # Parameters
    ratio::Int32,
    nfree::Int32,
    nboundary::Int32,
    numNodes::Int32,
    max_nnz_ii::Int32,
    max_nnz_b::Int32,
    use_varying_coeffs::Bool)

    g = CG.this_thread_block()
    tid = threadIdx().x
    blocksize = Int32(blockDim().x)
    iel = Int32(blockIdx().x)

    # ========================================================================
    # Shared Memory Allocation (only PCG workspace)
    # ========================================================================
    # Layout: diag_ii(nfree) | r(nfree) | z(nfree) | p(nfree) | Ap(nfree) | reduce(blocksize)

    shmem = CuDynamicSharedArray(Float64, 5 * nfree + blocksize)

    off_diag = Int32(0)
    off_r = nfree
    off_z = Int32(2) * nfree
    off_p = Int32(3) * nfree
    off_Ap = Int32(4) * nfree
    off_reduce = Int32(5) * nfree

    # ========================================================================
    # Load Coarse Element Corner Coordinates
    # ========================================================================

    # Use first 8 positions of shared memory temporarily
    if tid <= Int32(4)
        @inbounds node = cell_to_node[iel, tid]
        @inbounds shmem[2*tid - 1] = vertex_x[node]
        @inbounds shmem[2*tid] = vertex_y[node]
    end
    CG.sync(g)

    @inbounds x_corner1 = shmem[1]
    @inbounds y_corner1 = shmem[2]
    @inbounds x_corner2 = shmem[3]
    @inbounds y_corner2 = shmem[4]
    @inbounds x_corner3 = shmem[5]
    @inbounds y_corner3 = shmem[6]
    @inbounds x_corner4 = shmem[7]
    @inbounds y_corner4 = shmem[8]
    CG.sync(g)

    hx = (x_corner2 - x_corner1) / Float64(ratio)
    hy = (y_corner4 - y_corner1) / Float64(ratio)

    # ========================================================================
    # Phase 1: Initialize Per-Element Workspaces (Global Memory)
    # ========================================================================

    # Zero out per-element arrays
    i = tid
    while i <= max_nnz_ii
        @inbounds nzval_ii_all[iel, i] = 0.0
        i += blocksize
    end
    i = tid
    while i <= max_nnz_b
        @inbounds matValues_b_all[iel, i] = 0.0
        i += blocksize
    end
    i = tid
    while i <= numNodes
        @inbounds rhs_fine_all[iel, i] = 0.0
        for ir = Int32(1):Int32(4)
            @inbounds phi_all[iel, i, ir] = 0.0
        end
        i += blocksize
    end
    i = tid
    while i <= nfree
        for ir = Int32(1):Int32(3)
            @inbounds btmp_all[iel, i, ir] = 0.0
        end
        i += blocksize
    end
    CG.sync(g)

    # Initialize boundary phi values from pre-computed array
    i = tid
    while i <= nboundary
        @inbounds gi = boundaryToGlobal[i]
        for ir = Int32(1):Int32(4)
            @inbounds phi_all[iel, gi, ir] = phi_boundary[i, ir]
        end
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 2: Parallel Fine Element Assembly
    # ========================================================================

    nfine_elements = ratio * ratio
    ielem = tid
    while ielem <= nfine_elements
        ix = (ielem - Int32(1)) % ratio
        iy = (ielem - Int32(1)) ÷ ratio

        x1 = x_corner1 + Float64(ix) * hx
        y1 = y_corner1 + Float64(iy) * hy
        x2 = x1 + hx
        y2 = y1
        x3 = x2
        y3 = y2 + hy
        x4 = x1
        y4 = y3

        # Evaluate coefficients at each quadrature point inside assemble_fine_element_device
        Ke11, Ke12, Ke13, Ke14, Ke21, Ke22, Ke23, Ke24, Ke31, Ke32, Ke33, Ke34, Ke41, Ke42, Ke43, Ke44, fe1, fe2, fe3, fe4 =
            assemble_fine_element_device(x1, y1, x2, y2, x3, y3, x4, y4, use_varying_coeffs)

        node1 = ix + iy * (ratio + Int32(1)) + Int32(1)
        node2 = node1 + Int32(1)
        node3 = node2 + (ratio + Int32(1))
        node4 = node1 + (ratio + Int32(1))

        Ke = (Ke11, Ke12, Ke13, Ke14, Ke21, Ke22, Ke23, Ke24, Ke31, Ke32, Ke33, Ke34, Ke41, Ke42, Ke43, Ke44)
        fe = (fe1, fe2, fe3, fe4)
        nodes = (node1, node2, node3, node4)

        # Accumulate RHS (atomic to global memory)
        for i = 1:4
            @inbounds gi = nodes[i]
            @inbounds CUDA.@atomic rhs_fine_all[iel, gi] += fe[i]
        end

        # Assemble into sparse matrices
        for i = 1:4
            @inbounds iGlobal = nodes[i]
            @inbounds iFree = globalToFree[iGlobal]

            if iFree != Int32(-1)
                for j = 1:4
                    @inbounds jGlobal = nodes[j]
                    @inbounds jFree = globalToFree[jGlobal]
                    @inbounds k_val = Ke[(i-1)*4 + j]

                    if jFree != Int32(-1)
                        # Interior-interior coupling
                        k = find_matrix_position_device(colptr_ii, rowidx_ii, iFree, jFree)
                        if k != Int32(-1)
                            @inbounds CUDA.@atomic nzval_ii_all[iel, k] += k_val
                        end
                    else
                        # Interior-boundary coupling -> contributes to RHS
                        for ir = 1:3
                            @inbounds phi_val = phi_all[iel, jGlobal, ir]
                            @inbounds CUDA.@atomic btmp_all[iel, iFree, ir] -= k_val * phi_val
                        end
                    end
                end
            else
                # Boundary row -> K_b matrix
                @inbounds iBoundary = globalToBoundary[iGlobal]
                for j = 1:4
                    @inbounds jGlobal = nodes[j]
                    @inbounds k_val = Ke[(i-1)*4 + j]
                    k = find_matrix_position_device(matRowPtr_b, matColIdx_b, iBoundary, jGlobal)
                    if k != Int32(-1)
                        @inbounds CUDA.@atomic matValues_b_all[iel, k] += k_val
                    end
                end
            end
        end

        ielem += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 3: Extract Diagonal (to shared memory for PCG)
    # ========================================================================

    i = tid
    while i <= nfree
        diag_val = 0.0
        @inbounds for k = colptr_ii[i]:(colptr_ii[i + 1] - Int32(1))
            if rowidx_ii[k] == i
                @inbounds diag_val = nzval_ii_all[iel, k]
                break
            end
        end
        @inbounds shmem[off_diag + i] = diag_val
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 4: Initialize PCG Solution Vectors
    # ========================================================================

    i = tid
    while i <= nfree
        @inbounds gi = freeToGlobal[i]
        iy_fine = (gi - Int32(1)) ÷ (ratio + Int32(1))
        ix_fine = (gi - Int32(1)) % (ratio + Int32(1))
        xi = Float64(ix_fine) / Float64(ratio)
        eta = Float64(iy_fine) / Float64(ratio)

        # Bilinear initial guess
        N1 = (1.0 - xi) * (1.0 - eta)
        N2 = xi * (1.0 - eta)
        N3 = xi * eta

        @inbounds utmp_all[iel, i, 1] = N1
        @inbounds utmp_all[iel, i, 2] = N2
        @inbounds utmp_all[iel, i, 3] = N3
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 5: Solve for Interior Basis Functions with Block-Cooperative PCG
    # ========================================================================

    # Create views for shared memory arrays
    diag_ii_view = @inbounds view(shmem, (off_diag + 1):(off_diag + nfree))
    r_view = @inbounds view(shmem, (off_r + 1):(off_r + nfree))
    z_view = @inbounds view(shmem, (off_z + 1):(off_z + nfree))
    p_view = @inbounds view(shmem, (off_p + 1):(off_p + nfree))
    Ap_view = @inbounds view(shmem, (off_Ap + 1):(off_Ap + nfree))
    reduce_view = @inbounds view(shmem, (off_reduce + 1):(off_reduce + blocksize))

    for ir = Int32(1):Int32(3)
        # Get views into global memory for this RHS
        x_view = @inbounds view(utmp_all, iel, :, ir)
        b_view = @inbounds view(btmp_all, iel, :, ir)
        nzval_view = @inbounds view(nzval_ii_all, iel, :)

        block_pcg_jacobi_hybrid!(g, x_view, nfree,
                                  colptr_ii, rowidx_ii, nzval_view,
                                  b_view,
                                  diag_ii_view,
                                  r_view, z_view, p_view, Ap_view,
                                  reduce_view,
                                  blocksize,
                                  PCG_TOL, Int32(PCG_MAXITER))
    end
    CG.sync(g)

    # ========================================================================
    # Phase 6: Reconstruct Full Basis Functions
    # ========================================================================

    i = tid
    while i <= nfree
        @inbounds gi = freeToGlobal[i]
        sum_val = 0.0
        for ir = 1:3
            @inbounds val = utmp_all[iel, i, ir]
            @inbounds phi_all[iel, gi, ir] = val
            sum_val += val
        end
        @inbounds phi_all[iel, gi, 4] = 1.0 - sum_val
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 7: Compute Coarse Element Matrices
    # ========================================================================

    # Compute fe_coarse = phi^T * rhs_fine using block reduction
    for ir = Int32(1):Int32(4)
        local_sum = 0.0
        i = tid
        while i <= numNodes
            @inbounds local_sum += phi_all[iel, i, ir] * rhs_fine_all[iel, i]
            i += blocksize
        end
        dot_val = block_reduce_sum(g, reduce_view, local_sum, blocksize)
        if tid == Int32(1)
            @inbounds fe_out[iel, ir] = dot_val
        end
        CG.sync(g)
    end

    # Compute Ke_coarse = phi^T * K_b * phi
    # Each thread handles some boundary rows
    # First zero out Ke_out for this element
    if tid <= Int32(16)
        @inbounds Ke_out[iel, tid] = 0.0
    end
    CG.sync(g)

    i_bnd = tid
    while i_bnd <= nboundary
        @inbounds iGlobal = boundaryToGlobal[i_bnd]
        @inbounds rowBegin = matRowPtr_b[i_bnd]
        @inbounds rowEnd = matRowPtr_b[i_bnd + 1] - Int32(1)

        for k = rowBegin:rowEnd
            @inbounds jGlobal = matColIdx_b[k]
            @inbounds k_val = matValues_b_all[iel, k]

            for ir = 1:4
                @inbounds phi_i = phi_all[iel, iGlobal, ir]
                phi_i_k = phi_i * k_val
                for jr = 1:4
                    @inbounds phi_j = phi_all[iel, jGlobal, jr]
                    @inbounds CUDA.@atomic Ke_out[iel, (ir-1)*4 + jr] += phi_i_k * phi_j
                end
            end
        end
        i_bnd += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 8: Copy phi to output
    # ========================================================================

    i = tid
    while i <= numNodes
        for ir = 1:4
            @inbounds phi_out[iel, i, ir] = phi_all[iel, i, ir]
        end
        i += blocksize
    end

    return nothing
end

# ============================================================================
# GPU Coarse Global Assembly Kernel (from mfem_cuda.jl)
# ============================================================================

# Binary search helper for sparse matrix lookups
@inline function binary_search_sparse(indices, val, start_pos, end_pos)
    left = start_pos
    right = end_pos
    while left <= right
        mid = (left + right) ÷ 2
        @inbounds idx_val = indices[mid]
        if idx_val == val
            return mid
        elseif idx_val < val
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return -1  # Not found
end

# CUDA kernel to assemble coarse global matrix and RHS from element
# matrices. Uses coloring to avoid race conditions (all elements in same
# color processed in parallel)
function assemble_coarse_color_kernel!(
    d_mat_values, d_rhs, d_fe_coarse, d_Ke_coarse,
    d_cell_to_node, d_col_ptr, d_row_idx, d_elements, num_elements
)
    # One thread per element in this color
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if idx <= num_elements
        # Get element ID from color list
        iel = @inbounds d_elements[idx]

        # Get element nodes (Q1 element has 4 nodes)
        node1 = @inbounds d_cell_to_node[iel, 1]
        node2 = @inbounds d_cell_to_node[iel, 2]
        node3 = @inbounds d_cell_to_node[iel, 3]
        node4 = @inbounds d_cell_to_node[iel, 4]
        nodes = (node1, node2, node3, node4)

        # Scatter RHS: rhs[gi] += fe_coarse[iel, i]
        for i = 1:4
            gi = @inbounds nodes[i]
            fe_val = @inbounds d_fe_coarse[iel, i]
            @inbounds d_rhs[gi] += fe_val
        end

        # Scatter stiffness: mat_values[k] += Ke_coarse[iel, i, j]
        # Loop over columns first to reuse col_ptr loads (4 loads instead of 16)
        for j = 1:4
            gj = @inbounds nodes[j]

            # Load column bounds once per column (not once per entry)
            col_start = @inbounds d_col_ptr[gj]
            col_end = @inbounds d_col_ptr[gj + 1] - 1

            for i = 1:4
                gi = @inbounds nodes[i]

                # Binary search for row gi in column gj
                k = binary_search_sparse(d_row_idx, gi, col_start, col_end)

                if k != -1
                    # Ke_coarse is stored row-major: (i-1)*4 + j
                    ke_idx = (i - 1) * 4 + j
                    ke_val = @inbounds d_Ke_coarse[iel, ke_idx]
                    @inbounds d_mat_values[k] += ke_val
                end
            end
        end
    end

    return nothing
end

# ============================================================================
# Pre-computation Functions (Host Side)
# ============================================================================

"""
Build DOF mappings for a given ratio.
These are the same for all coarse elements.
"""
function build_dof_mappings(ratio::Int)
    numNodes = (ratio + 1)^2
    nfree = (ratio - 1)^2
    nboundary = 4 * ratio

    globalToFree = fill(Int32(-1), numNodes)
    freeToGlobal = zeros(Int32, nfree)
    globalToBoundary = fill(Int32(-1), numNodes)
    boundaryToGlobal = zeros(Int32, nboundary)

    nfree_count = 0
    nboundary_count = 0

    for iy = 0:ratio
        for ix = 0:ratio
            nodeID = ix + iy * (ratio + 1) + 1
            if ix > 0 && ix < ratio && iy > 0 && iy < ratio
                nfree_count += 1
                freeToGlobal[nfree_count] = Int32(nodeID)
                globalToFree[nodeID] = Int32(nfree_count)
            else
                nboundary_count += 1
                boundaryToGlobal[nboundary_count] = Int32(nodeID)
                globalToBoundary[nodeID] = Int32(nboundary_count)
            end
        end
    end

    return globalToFree, freeToGlobal, globalToBoundary, boundaryToGlobal
end

"""
Build K_ii sparsity pattern (CSC format).
Same for all coarse elements with given ratio.
"""
function build_kii_sparsity(ratio::Int, globalToFree::Vector{Int32}, freeToGlobal::Vector{Int32})
    nfree = length(freeToGlobal)

    # Count non-zeros per column
    colptr = zeros(Int32, nfree + 1)
    colptr[1] = Int32(1)

    for i_free = 1:nfree
        iGlobal = freeToGlobal[i_free]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)

        count = 1  # diagonal
        if iy > 1
            if ix > 1; count += 1; end
            count += 1
            if ix < ratio - 1; count += 1; end
        end
        if ix > 1; count += 1; end
        if ix < ratio - 1; count += 1; end
        if iy < ratio - 1
            if ix > 1; count += 1; end
            count += 1
            if ix < ratio - 1; count += 1; end
        end
        colptr[i_free + 1] = colptr[i_free] + Int32(count)
    end

    nnz = colptr[end] - 1
    rowidx = zeros(Int32, nnz)

    # Fill row indices
    for i_free = 1:nfree
        iGlobal = freeToGlobal[i_free]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)
        offset = colptr[i_free]

        # South neighbors
        if iy > 1
            if ix > 1
                jGlobal = iGlobal - 1 - (ratio + 1)
                jFree = globalToFree[jGlobal]
                if jFree != -1
                    rowidx[offset] = jFree; offset += 1
                end
            end
            jGlobal = iGlobal - (ratio + 1)
            jFree = globalToFree[jGlobal]
            if jFree != -1
                rowidx[offset] = jFree; offset += 1
            end
            if ix < ratio - 1
                jGlobal = iGlobal + 1 - (ratio + 1)
                jFree = globalToFree[jGlobal]
                if jFree != -1
                    rowidx[offset] = jFree; offset += 1
                end
            end
        end
        # West neighbor
        if ix > 1
            jGlobal = iGlobal - 1
            jFree = globalToFree[jGlobal]
            if jFree != -1
                rowidx[offset] = jFree; offset += 1
            end
        end
        # Diagonal
        rowidx[offset] = Int32(i_free); offset += 1
        # East neighbor
        if ix < ratio - 1
            jGlobal = iGlobal + 1
            jFree = globalToFree[jGlobal]
            if jFree != -1
                rowidx[offset] = jFree; offset += 1
            end
        end
        # North neighbors
        if iy < ratio - 1
            if ix > 1
                jGlobal = iGlobal - 1 + (ratio + 1)
                jFree = globalToFree[jGlobal]
                if jFree != -1
                    rowidx[offset] = jFree; offset += 1
                end
            end
            jGlobal = iGlobal + (ratio + 1)
            jFree = globalToFree[jGlobal]
            if jFree != -1
                rowidx[offset] = jFree; offset += 1
            end
            if ix < ratio - 1
                jGlobal = iGlobal + 1 + (ratio + 1)
                jFree = globalToFree[jGlobal]
                if jFree != -1
                    rowidx[offset] = jFree; offset += 1
                end
            end
        end
    end

    return colptr, rowidx
end

"""
Build K_b sparsity pattern (CSR format, rows indexed by boundary DOF, cols by global node).
"""
function build_kb_sparsity(ratio::Int, boundaryToGlobal::Vector{Int32})
    nboundary = length(boundaryToGlobal)
    numNodes = (ratio + 1)^2

    rowptr = zeros(Int32, nboundary + 1)
    rowptr[1] = Int32(1)

    for i_bnd = 1:nboundary
        iGlobal = boundaryToGlobal[i_bnd]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)

        count = 1  # self
        if ix > 0; count += 1; end
        if ix < ratio; count += 1; end
        if iy > 0
            count += 1
            if ix > 0; count += 1; end
            if ix < ratio; count += 1; end
        end
        if iy < ratio
            count += 1
            if ix > 0; count += 1; end
            if ix < ratio; count += 1; end
        end
        rowptr[i_bnd + 1] = rowptr[i_bnd] + Int32(count)
    end

    nnz = rowptr[end] - 1
    colidx = zeros(Int32, nnz)

    for i_bnd = 1:nboundary
        iGlobal = boundaryToGlobal[i_bnd]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)
        offset = rowptr[i_bnd]

        if iy > 0
            if ix > 0
                colidx[offset] = Int32(iGlobal - 1 - (ratio + 1)); offset += 1
            end
            colidx[offset] = Int32(iGlobal - (ratio + 1)); offset += 1
            if ix < ratio
                colidx[offset] = Int32(iGlobal + 1 - (ratio + 1)); offset += 1
            end
        end
        if ix > 0
            colidx[offset] = Int32(iGlobal - 1); offset += 1
        end
        colidx[offset] = Int32(iGlobal); offset += 1
        if ix < ratio
            colidx[offset] = Int32(iGlobal + 1); offset += 1
        end
        if iy < ratio
            if ix > 0
                colidx[offset] = Int32(iGlobal - 1 + (ratio + 1)); offset += 1
            end
            colidx[offset] = Int32(iGlobal + (ratio + 1)); offset += 1
            if ix < ratio
                colidx[offset] = Int32(iGlobal + 1 + (ratio + 1)); offset += 1
            end
        end
    end

    return rowptr, colidx
end

"""
Build boundary phi values (same for all elements with given ratio).
"""
function build_boundary_phi(ratio::Int, boundaryToGlobal::Vector{Int32})
    nboundary = length(boundaryToGlobal)
    numNodes = (ratio + 1)^2

    phi_boundary = zeros(Float64, nboundary, 4)

    # Set corner values
    corner1 = 1
    corner2 = ratio + 1
    corner3 = numNodes
    corner4 = (ratio + 1) * ratio + 1

    # Find boundary indices for corners
    for i = 1:nboundary
        gi = boundaryToGlobal[i]
        if gi == corner1
            phi_boundary[i, 1] = 1.0
        elseif gi == corner2
            phi_boundary[i, 2] = 1.0
        elseif gi == corner3
            phi_boundary[i, 3] = 1.0
        elseif gi == corner4
            phi_boundary[i, 4] = 1.0
        end
    end

    # Interpolate along edges
    for i = 1:nboundary
        gi = boundaryToGlobal[i]
        ix = (gi - 1) % (ratio + 1)
        iy = (gi - 1) ÷ (ratio + 1)

        # Skip corners (already set)
        if (ix == 0 && iy == 0) || (ix == ratio && iy == 0) ||
           (ix == ratio && iy == ratio) || (ix == 0 && iy == ratio)
            continue
        end

        if ix == 0  # Left edge
            s = Float64(iy) / Float64(ratio)
            phi_boundary[i, 1] = 1.0 - s
            phi_boundary[i, 4] = s
        elseif ix == ratio  # Right edge
            s = Float64(iy) / Float64(ratio)
            phi_boundary[i, 2] = 1.0 - s
            phi_boundary[i, 3] = s
        elseif iy == 0  # Bottom edge
            s = Float64(ix) / Float64(ratio)
            phi_boundary[i, 1] = 1.0 - s
            phi_boundary[i, 2] = s
        elseif iy == ratio  # Top edge
            s = Float64(ix) / Float64(ratio)
            phi_boundary[i, 4] = 1.0 - s
            phi_boundary[i, 3] = s
        end
    end

    return phi_boundary
end

# ============================================================================
# Host-Side Interface
# ============================================================================

struct CudaMesh
    vertex_x::CuVector{Float64}
    vertex_y::CuVector{Float64}
    cell_to_node::CuMatrix{Int32}
    nel::Int
    nnodes::Int
end

function CudaMesh(vertex_x::Vector{Float64}, vertex_y::Vector{Float64},
                  cell_to_node::Matrix{Int})
    nel = size(cell_to_node, 1)
    nnodes = length(vertex_x)
    CudaMesh(
        CuVector(vertex_x),
        CuVector(vertex_y),
        CuMatrix(Int32.(cell_to_node)),
        nel,
        nnodes
    )
end

"""
Pre-computed data for a given ratio (shared across all elements).
"""
struct MFEMPrecomputed
    ratio::Int
    nfree::Int
    nboundary::Int
    numNodes::Int
    max_nnz_ii::Int
    max_nnz_b::Int
    # DOF mappings (on GPU)
    globalToFree::CuVector{Int32}
    freeToGlobal::CuVector{Int32}
    globalToBoundary::CuVector{Int32}
    boundaryToGlobal::CuVector{Int32}
    # Sparsity patterns (on GPU)
    colptr_ii::CuVector{Int32}
    rowidx_ii::CuVector{Int32}
    matRowPtr_b::CuVector{Int32}
    matColIdx_b::CuVector{Int32}
    # Boundary phi values (on GPU)
    phi_boundary::CuMatrix{Float64}
end

function MFEMPrecomputed(ratio::Int)
    numNodes = (ratio + 1)^2
    nfree = (ratio - 1)^2
    nboundary = 4 * ratio

    # Build on CPU
    globalToFree, freeToGlobal, globalToBoundary, boundaryToGlobal = build_dof_mappings(ratio)
    colptr_ii, rowidx_ii = build_kii_sparsity(ratio, globalToFree, freeToGlobal)
    matRowPtr_b, matColIdx_b = build_kb_sparsity(ratio, boundaryToGlobal)
    phi_boundary = build_boundary_phi(ratio, boundaryToGlobal)

    max_nnz_ii = length(rowidx_ii)
    max_nnz_b = length(matColIdx_b)

    MFEMPrecomputed(
        ratio, nfree, nboundary, numNodes, max_nnz_ii, max_nnz_b,
        CuVector(globalToFree),
        CuVector(freeToGlobal),
        CuVector(globalToBoundary),
        CuVector(boundaryToGlobal),
        CuVector(colptr_ii),
        CuVector(rowidx_ii),
        CuVector(matRowPtr_b),
        CuVector(matColIdx_b),
        CuMatrix(phi_boundary)
    )
end

function calculate_shmem_size(nfree::Int, blocksize::Int)
    # Shared memory: diag(nfree) + r(nfree) + z(nfree) + p(nfree) + Ap(nfree) + reduce(blocksize)
    return (5 * nfree + blocksize) * sizeof(Float64)
end

function assemble_mfem_cuda_block(mesh::CudaMesh, precomp::MFEMPrecomputed;
                                   blocksize::Int=256,
                                   use_varying_coeffs::Bool=true)
    nel = mesh.nel
    ratio = precomp.ratio
    nfree = precomp.nfree
    nboundary = precomp.nboundary
    numNodes = precomp.numNodes
    max_nnz_ii = precomp.max_nnz_ii
    max_nnz_b = precomp.max_nnz_b

    # Allocate outputs
    Ke_out = CUDA.zeros(Float64, nel, 16)
    fe_out = CUDA.zeros(Float64, nel, 4)
    phi_out = CUDA.zeros(Float64, nel, numNodes, 4)

    # Allocate per-element workspaces (global memory)
    nzval_ii_all = CUDA.zeros(Float64, nel, max_nnz_ii)
    matValues_b_all = CUDA.zeros(Float64, nel, max_nnz_b)
    phi_all = CUDA.zeros(Float64, nel, numNodes, 4)
    rhs_fine_all = CUDA.zeros(Float64, nel, numNodes)
    btmp_all = CUDA.zeros(Float64, nel, nfree, 3)
    utmp_all = CUDA.zeros(Float64, nel, nfree, 3)

    shmem_size = calculate_shmem_size(nfree, blocksize)

    @cuda threads=blocksize blocks=nel shmem=shmem_size mfem_element_kernel_hybrid!(
        Ke_out, fe_out, phi_out,
        mesh.vertex_x, mesh.vertex_y, mesh.cell_to_node,
        precomp.colptr_ii, precomp.rowidx_ii,
        precomp.matRowPtr_b, precomp.matColIdx_b,
        precomp.globalToFree, precomp.freeToGlobal,
        precomp.globalToBoundary, precomp.boundaryToGlobal,
        precomp.phi_boundary,
        nzval_ii_all, matValues_b_all, phi_all,
        rhs_fine_all, btmp_all, utmp_all,
        Int32(ratio), Int32(nfree), Int32(nboundary), Int32(numNodes),
        Int32(max_nnz_ii), Int32(max_nnz_b),
        use_varying_coeffs
    )

    return Ke_out, fe_out, phi_out
end

# Convenience function that creates precomputed data
function assemble_mfem_cuda_block(mesh::CudaMesh, ratio::Int; kwargs...)
    precomp = MFEMPrecomputed(ratio)
    return assemble_mfem_cuda_block(mesh, precomp; kwargs...)
end

# ============================================================================
# Test Functions
# ============================================================================

function generate_test_mesh(nx::Int, ny::Int)
    nv = (nx + 1) * (ny + 1)
    nel = nx * ny

    vertex_x = zeros(nv)
    vertex_y = zeros(nv)
    cell_to_node = zeros(Int, nel, 4)

    for j = 0:ny, i = 0:nx
        idx = i + 1 + j * (nx + 1)
        vertex_x[idx] = Float64(i) / nx
        vertex_y[idx] = Float64(j) / ny
    end

    for iel = 1:nel
        i = (iel - 1) % nx
        j = (iel - 1) ÷ nx
        n1 = j * (nx + 1) + i + 1
        cell_to_node[iel, 1] = n1
        cell_to_node[iel, 2] = n1 + 1
        cell_to_node[iel, 3] = n1 + 1 + (nx + 1)
        cell_to_node[iel, 4] = n1 + (nx + 1)
    end

    return vertex_x, vertex_y, cell_to_node
end

function test_mfem_cuda_block()
    println("="^70)
    println("CUDA MFEM Assembly Test (Hybrid Memory Version)")
    println("="^70)
    println()

    if !CUDA.functional()
        println("CUDA not available, skipping test")
        return
    end

    println("CUDA Device: ", CUDA.name(CUDA.device()))
    println()

    nx, ny = 4, 4
    ratio = 8

    println("Mesh: $nx × $ny coarse elements")
    println("Ratio: $ratio ($(ratio)×$(ratio) fine elements per coarse)")
    println("Effective fine mesh: $(nx*ratio) × $(ny*ratio)")
    println()

    # Pre-compute sparsity patterns
    println("Pre-computing sparsity patterns...")
    precomp = MFEMPrecomputed(ratio)
    println("  nfree = $(precomp.nfree)")
    println("  nboundary = $(precomp.nboundary)")
    println("  nnz(K_ii) = $(precomp.max_nnz_ii)")
    println("  nnz(K_b) = $(precomp.max_nnz_b)")
    println()

    vertex_x, vertex_y, cell_to_node = generate_test_mesh(nx, ny)
    mesh = CudaMesh(vertex_x, vertex_y, cell_to_node)

    shmem_size = calculate_shmem_size(precomp.nfree, 256)
    println("Shared memory per block: $(shmem_size) bytes ($(shmem_size ÷ 1024) KB)")
    println()

    println("Running warm-up...")
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp)
    CUDA.synchronize()

    println("Running timed assembly...")
    t_start = time()
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp)
    CUDA.synchronize()
    t_elapsed = time() - t_start

    println()
    println("Assembly time: $(t_elapsed * 1000) ms")
    println()

    Ke_cpu = Array(Ke_out)
    fe_cpu = Array(fe_out)

    println("Sample results (element 1):")
    println("  Ke[1,1:4] = ", Ke_cpu[1, 1:4])
    println("  fe[1,:] = ", fe_cpu[1, :])
    println()

    Ke1 = reshape(Ke_cpu[1, :], 4, 4)
    sym_error = maximum(abs.(Ke1 - Ke1'))
    println("Ke symmetry error: $sym_error")

    eigvals_Ke = eigvals(Ke1)
    println("Ke eigenvalues: ", eigvals_Ke)
    println()

    println("Test completed successfully!")
end

function benchmark_mfem_cuda_block(nx::Int=16, ny::Int=16, ratio::Int=8; nruns::Int=10)
    println("="^70)
    println("CUDA MFEM Assembly Benchmark (Hybrid Memory Version)")
    println("="^70)
    println()

    if !CUDA.functional()
        println("CUDA not available, skipping benchmark")
        return
    end

    println("CUDA Device: ", CUDA.name(CUDA.device()))
    println()

    println("Mesh: $nx × $ny coarse elements")
    println("Ratio: $ratio")
    println("Effective fine mesh: $(nx*ratio) × $(ny*ratio)")
    println("Number of coarse elements: $(nx*ny)")
    println()

    precomp = MFEMPrecomputed(ratio)
    println("Pre-computed data:")
    println("  nfree = $(precomp.nfree)")
    println("  nnz(K_ii) = $(precomp.max_nnz_ii)")
    println()

    vertex_x, vertex_y, cell_to_node = generate_test_mesh(nx, ny)
    mesh = CudaMesh(vertex_x, vertex_y, cell_to_node)

    # Warm up
    for _ = 1:3
        Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp)
        CUDA.synchronize()
    end

    times = Float64[]
    for _ = 1:nruns
        t = CUDA.@elapsed begin
            Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp)
            CUDA.synchronize()
        end
        push!(times, t)
    end

    mean_time = sum(times) / length(times)
    min_time = minimum(times)
    max_time = maximum(times)

    println("Results over $nruns runs:")
    println("  Mean time: $(mean_time * 1000) ms")
    println("  Min time:  $(min_time * 1000) ms")
    println("  Max time:  $(max_time * 1000) ms")
    println("  Elements/sec: $(nx*ny / mean_time)")
    println()
end

function test_large_ratio()
    println("="^70)
    println("Testing Large Ratio Support (Hybrid Memory)")
    println("="^70)
    println()

    if !CUDA.functional()
        println("CUDA not available, skipping test")
        return
    end

    for ratio in [8, 16, 32, 64]
        println("Testing ratio = $ratio")
        precomp = MFEMPrecomputed(ratio)
        shmem_size = calculate_shmem_size(precomp.nfree, 256)
        println("  nfree = $(precomp.nfree)")
        println("  Shared memory: $(shmem_size ÷ 1024) KB")

        nx, ny = 4, 4
        vertex_x, vertex_y, cell_to_node = generate_test_mesh(nx, ny)
        mesh = CudaMesh(vertex_x, vertex_y, cell_to_node)

        try
            Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp)
            CUDA.synchronize()

            Ke_cpu = Array(Ke_out)
            Ke1 = reshape(Ke_cpu[1, :], 4, 4)
            sym_error = maximum(abs.(Ke1 - Ke1'))
            println("  Symmetry error: $sym_error")
            println("  ✓ Success")
        catch e
            println("  ✗ Failed: $e")
        end
        println()
    end
end

# ============================================================================
# Global Assembly and Solve (CPU)
# ============================================================================

"""
Build node-to-node sparsity pattern for a structured Q1 mesh.
Returns row_ptr, col_idx (CSR format) for the coarse mesh.
"""
function build_n2n_sparsity(nx::Int, ny::Int)
    nnodes = (nx + 1) * (ny + 1)

    # Count neighbors for each node
    row_ptr = zeros(Int, nnodes + 1)
    row_ptr[1] = 1

    for j = 0:ny
        for i = 0:nx
            node = i + 1 + j * (nx + 1)
            count = 1  # self
            if i > 0; count += 1; end
            if i < nx; count += 1; end
            if j > 0; count += 1; end
            if j < ny; count += 1; end
            if i > 0 && j > 0; count += 1; end
            if i < nx && j > 0; count += 1; end
            if i > 0 && j < ny; count += 1; end
            if i < nx && j < ny; count += 1; end
            row_ptr[node + 1] = row_ptr[node] + count
        end
    end

    nnz = row_ptr[end] - 1
    col_idx = zeros(Int, nnz)

    # Fill column indices (sorted order)
    for j = 0:ny
        for i = 0:nx
            node = i + 1 + j * (nx + 1)
            offset = row_ptr[node]

            # Add neighbors in sorted order
            neighbors = Int[]
            if j > 0 && i > 0; push!(neighbors, node - (nx + 1) - 1); end
            if j > 0; push!(neighbors, node - (nx + 1)); end
            if j > 0 && i < nx; push!(neighbors, node - (nx + 1) + 1); end
            if i > 0; push!(neighbors, node - 1); end
            push!(neighbors, node)  # self
            if i < nx; push!(neighbors, node + 1); end
            if j < ny && i > 0; push!(neighbors, node + (nx + 1) - 1); end
            if j < ny; push!(neighbors, node + (nx + 1)); end
            if j < ny && i < nx; push!(neighbors, node + (nx + 1) + 1); end

            sort!(neighbors)
            for (k, n) in enumerate(neighbors)
                col_idx[offset + k - 1] = n
            end
        end
    end

    return row_ptr, col_idx
end

"""
Get boundary nodes for a structured mesh on [0,1]^2.
"""
function get_boundary_nodes(nx::Int, ny::Int)
    boundary = Int[]

    # Bottom edge (j=0)
    for i = 0:nx
        push!(boundary, i + 1)
    end
    # Top edge (j=ny)
    for i = 0:nx
        push!(boundary, i + 1 + ny * (nx + 1))
    end
    # Left edge (i=0), excluding corners
    for j = 1:ny-1
        push!(boundary, 1 + j * (nx + 1))
    end
    # Right edge (i=nx), excluding corners
    for j = 1:ny-1
        push!(boundary, nx + 1 + j * (nx + 1))
    end

    return sort!(unique!(boundary))
end

"""
Find position in CSR matrix using binary search.
"""
function find_matrix_position(row_ptr::Vector{Int}, col_idx::Vector{Int}, row::Int, col::Int)
    left = row_ptr[row]
    right = row_ptr[row + 1] - 1

    while left <= right
        mid = (left + right) >> 1
        if col_idx[mid] == col
            return mid
        elseif col_idx[mid] < col
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return -1
end

"""
Scatter element matrices to global sparse matrix.
"""
function scatter_to_global!(mat_values::Vector{Float64}, rhs::Vector{Float64},
                            Ke_cpu::Matrix{Float64}, fe_cpu::Matrix{Float64},
                            cell_to_node::Matrix{Int},
                            row_ptr::Vector{Int}, col_idx::Vector{Int})
    nel = size(cell_to_node, 1)

    for iel = 1:nel
        for i = 1:4
            gi = cell_to_node[iel, i]
            rhs[gi] += fe_cpu[iel, i]
            for j = 1:4
                gj = cell_to_node[iel, j]
                k = find_matrix_position(row_ptr, col_idx, gi, gj)
                if k != -1
                    mat_values[k] += Ke_cpu[iel, (i-1)*4 + j]
                end
            end
        end
    end
end

"""
Apply homogeneous Dirichlet boundary conditions.
Returns reduced system.
"""
function apply_boundary_conditions(row_ptr::Vector{Int}, col_idx::Vector{Int},
                                   mat_values::Vector{Float64}, rhs::Vector{Float64},
                                   boundary_nodes::Vector{Int})
    nnodes = length(rhs)

    # Build DOF mapping
    global_to_free = fill(-1, nnodes)
    free_to_global = Int[]

    bdy_set = Set(boundary_nodes)
    for i = 1:nnodes
        if !(i in bdy_set)
            push!(free_to_global, i)
            global_to_free[i] = length(free_to_global)
        end
    end

    nfree = length(free_to_global)

    # Build reduced system
    reduced_row_ptr = zeros(Int, nfree + 1)
    reduced_col_idx = Int[]
    reduced_values = Float64[]
    reduced_rhs = zeros(nfree)

    reduced_row_ptr[1] = 1
    for i = 1:nfree
        gi = free_to_global[i]
        reduced_rhs[i] = rhs[gi]
        for k = row_ptr[gi]:row_ptr[gi + 1] - 1
            gj = col_idx[k]
            if global_to_free[gj] != -1
                push!(reduced_col_idx, global_to_free[gj])
                push!(reduced_values, mat_values[k])
            end
        end
        reduced_row_ptr[i + 1] = length(reduced_col_idx) + 1
    end

    return reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global
end

"""
Reconstruct fine-scale solution using MFEM basis functions.
"""
function reconstruct_fine_solution(nx::Int, ny::Int, ratio::Int,
                                   coarse_solution::Vector{Float64},
                                   phi_cpu::Array{Float64, 3},
                                   cell_to_node::Matrix{Int})
    nel = nx * ny
    nx_fine = nx * ratio
    ny_fine = ny * ratio
    nfine = (nx_fine + 1) * (ny_fine + 1)
    numNodes_local = (ratio + 1)^2

    fine_x = zeros(nfine)
    fine_y = zeros(nfine)
    fine_u = zeros(nfine)

    # Build fine mesh coordinates
    for j = 0:ny_fine
        for i = 0:nx_fine
            idx = (i + 1) + j * (nx_fine + 1)
            fine_x[idx] = Float64(i) / nx_fine
            fine_y[idx] = Float64(j) / ny_fine
        end
    end

    # Reconstruct solution element by element
    for iel = 1:nel
        # Coarse element position
        ix_coarse = (iel - 1) % nx
        iy_coarse = (iel - 1) ÷ nx

        # Coarse solution at element corners
        u_coarse = zeros(4)
        for i = 1:4
            gi = cell_to_node[iel, i]
            u_coarse[i] = coarse_solution[gi]
        end

        # Loop over fine nodes in this coarse element
        for iy_local = 0:ratio
            for ix_local = 0:ratio
                # Local node index
                local_node = ix_local + iy_local * (ratio + 1) + 1

                # Global fine node index
                ix_fine = ix_coarse * ratio + ix_local
                iy_fine = iy_coarse * ratio + iy_local
                fine_node = (ix_fine + 1) + iy_fine * (nx_fine + 1)

                # Interpolate using MFEM basis functions
                u_fine = 0.0
                for k = 1:4
                    u_fine += phi_cpu[iel, local_node, k] * u_coarse[k]
                end
                fine_u[fine_node] = u_fine
            end
        end
    end

    return fine_x, fine_y, fine_u
end

"""
Write solution to file.
"""
function write_solution(filename::String, vertex_x::Vector{Float64},
                        vertex_y::Vector{Float64}, solution::Vector{Float64})
    open(filename, "w") do io
        println(io, "# x y u")
        for i = 1:length(solution)
            @printf(io, "%.6f %.6f %.6e\n", vertex_x[i], vertex_y[i], solution[i])
        end
    end
end

# ============================================================================
# Main Function (CLI Interface)
# ============================================================================

function main()
    println("="^70)
    println("Julia 2D MFEM Assembly (CUDA - One Block Per Coarse Element)")
    println("="^70)
    println()
    println("Usage: julia mfem_cuda_block.jl [nx] [ny] [ratio]")
    println("  nx:     Coarse mesh elements in x-direction (default: 48)")
    println("  ny:     Coarse mesh elements in y-direction (default: nx)")
    println("  ratio:  Fine grid refinement ratio (default: 64)")
    println()

    if !CUDA.functional()
        println("ERROR: CUDA not available")
        return
    end

    # Device info (similar to thread info in CPU version)
    println("CUDA Device: ", CUDA.name(CUDA.device()))
    println()

    # ========================================================================
    # Problem setup
    # ========================================================================

    USE_VARYING_COEFFICIENTS = true

    # ========================================================================
    # Parse command-line arguments
    # ========================================================================

    nx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 48
    ny = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nx
    ratio = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 64

    println("Generating coarse mesh: $nx x $ny MFEM elements")
    println("MFEM ratio: $ratio (each coarse element has $(ratio)x$(ratio) fine elements)")
    println("Effective fine resolution: $(nx*ratio) x $(ny*ratio)")
    if USE_VARYING_COEFFICIENTS
        println("Diffusion coefficients: VARYING (1 + 100 * cos(150*x)^2 * sin(150*y)^2 terms)")
    else
        println("Diffusion coefficients: CONSTANT (ax=ay=1)")
    end
    println()

    # ========================================================================
    # Mesh generation
    # ========================================================================

    println("Generating mesh...")
    t0 = time()
    # Use FEMBase mesh generation to get Mesh struct with boundary nodes
    fem_mesh = FEMBase.generate_mesh(nx, ny)
    mesh = CudaMesh(fem_mesh.vertex_x, fem_mesh.vertex_y, fem_mesh.cell_to_node)
    mesh_time = time() - t0

    println("  Number of coarse elements: ", mesh.nel)
    println("  Number of coarse nodes:    ", mesh.nnodes)
    println("  Mesh generation:           ", @sprintf("%.2f ms", mesh_time * 1000))
    println()

    # ========================================================================
    # Build connectivity and coloring (for GPU global assembly)
    # ========================================================================

    println("Building mesh connectivity and coloring...")
    t0 = time()
    n2n_row_ptr, n2n_col_idx, e2e_colors = FEMBase.build_mesh_connectivity(fem_mesh)
    conn_time = time() - t0

    nnz_coarse = length(n2n_col_idx)
    println("  Coarse matrix size: $(mesh.nnodes) x $(mesh.nnodes)")
    println("  Non-zeros:          $nnz_coarse")
    println("  Number of colors:   $(length(e2e_colors))")
    println("  Connectivity time:  ", @sprintf("%.2f ms", conn_time * 1000))
    println()

    # ========================================================================
    # Pre-compute sparsity patterns (analogous to connectivity/coloring in CPU version)
    # ========================================================================

    println("Pre-computing element sparsity patterns...")
    CUDA.synchronize()
    t0 = time()
    precomp = MFEMPrecomputed(ratio)
    CUDA.synchronize()
    precomp_time = time() - t0

    shmem_size = calculate_shmem_size(precomp.nfree, 256)
    println("  Interior DOFs per element: ", precomp.nfree)
    println("  Boundary DOFs per element: ", precomp.nboundary)
    println("  nnz(K_ii): ", precomp.max_nnz_ii, ", nnz(K_b): ", precomp.max_nnz_b)
    println("  Shared memory per block:   $(shmem_size) bytes ($(shmem_size ÷ 1024) KB)")
    println("  Pre-computation time:      ", @sprintf("%.2f ms", precomp_time * 1000))
    println()

    # ========================================================================
    # Assembly
    # ========================================================================

    println("="^70)
    println("Starting MFEM Assembly")
    println("="^70)

    # Warm-up run to eliminate JIT compilation overhead
    println("Running JIT warm-up...")
    t0 = time()
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp; use_varying_coeffs=USE_VARYING_COEFFICIENTS)
    CUDA.synchronize()
    warmup_time = time() - t0
    println("  Warm-up time: ", @sprintf("%.2f ms", warmup_time * 1000))
    println()

    # Timed run
    t0 = time()
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, precomp; use_varying_coeffs=USE_VARYING_COEFFICIENTS)
    CUDA.synchronize()
    assembly_time = time() - t0

    # ========================================================================
    # GPU Assembly for the coarse global matrix (like mfem_cuda.jl)
    # ========================================================================

    println("Assembling coarse global matrix on GPU using coloring...")

    # Allocate global matrix and RHS on GPU
    nnodes = mesh.nnodes
    nnz = length(n2n_col_idx)
    d_mat_values = CUDA.zeros(Float64, nnz)
    d_rhs = CUDA.zeros(Float64, nnodes)

    # Upload sparse matrix structure to GPU (CSC format: col_ptr, row_idx)
    d_col_ptr = cu(n2n_row_ptr)  # Note: n2n_row_ptr is col_ptr in CSC
    d_row_idx = cu(n2n_col_idx)  # Note: n2n_col_idx is row_idx in CSC
    d_cell_to_node = cu(Int32.(fem_mesh.cell_to_node))

    # Pre-upload all color element lists (avoid transfers in timed section)
    d_color_elements = [cu(Int32.(elements)) for elements in e2e_colors]

    # Warm-up: launch kernel once to trigger JIT compilation before timing
    if !isempty(e2e_colors)
        d_warmup = d_color_elements[1]
        num_warmup = length(e2e_colors[1])
        threads_per_block = 256
        num_blocks = cld(num_warmup, threads_per_block)
        @cuda threads=threads_per_block blocks=num_blocks (
            assemble_coarse_color_kernel!(
                d_mat_values, d_rhs, fe_out, Ke_out,
                d_cell_to_node, d_col_ptr, d_row_idx,
                d_warmup, num_warmup
            )
        )
        CUDA.synchronize()
        # Reset matrix values after warm-up (they got modified)
        CUDA.fill!(d_mat_values, 0.0)
        CUDA.fill!(d_rhs, 0.0)
    end

    # Assemble by color to avoid race conditions
    # Color element lists already pre-uploaded, no transfers in timed section
    CUDA.synchronize()
    t0 = time()
    for (ic, d_elements) in enumerate(d_color_elements)
        num_elements = length(e2e_colors[ic])

        # Launch kernel: one thread per element in this color
        threads_per_block = 256
        num_blocks = cld(num_elements, threads_per_block)

        @cuda threads=threads_per_block blocks=num_blocks (
            assemble_coarse_color_kernel!(
                d_mat_values, d_rhs, fe_out, Ke_out,
                d_cell_to_node, d_col_ptr, d_row_idx,
                d_elements, num_elements
            )
        )
    end
    CUDA.synchronize()
    gpu_scatter_time = time() - t0

    println("  Number of colors: $(length(e2e_colors))")
    println("  GPU scatter time: ", @sprintf("%.2f ms", gpu_scatter_time * 1000))

    # Transfer assembled matrix and RHS to CPU (like mfem_cuda.jl)
    CUDA.synchronize()
    t0 = time()
    mat_values = Array(d_mat_values)
    rhs = Array(d_rhs)
    CUDA.synchronize()
    transfer_time = time() - t0

    phi_cpu = Array(phi_out)

    println("MFEM assembly complete")
    println("  CUDA kernel time:  ", @sprintf("%.2f ms", assembly_time * 1000))
    println("  GPU scatter time:  ", @sprintf("%.2f ms", gpu_scatter_time * 1000))
    println("  Matrix size: $nnodes x $nnodes, nnz: $nnz")
    println("Transfer time =     ", @sprintf("%.2f ms", transfer_time * 1000))
    println()

    # ========================================================================
    # Apply boundary conditions
    # ========================================================================

    println("Applying boundary conditions...")
    t0 = time()
    # Use boundary nodes from FEMBase mesh
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global =
        FEMBase.apply_boundary_conditions(n2n_row_ptr, n2n_col_idx, mat_values, rhs, fem_mesh.boundary_nodes)
    bc_time = time() - t0

    nfree = length(free_to_global)
    println("  Total DOFs:    $nnodes")
    println("  Boundary DOFs: $(length(fem_mesh.boundary_nodes))")
    println("  Free DOFs:     $nfree")
    println("  Reduced nnz:   $(length(reduced_col_idx))")
    println("  BC time:       ", @sprintf("%.2f ms", bc_time * 1000))
    println()

    # ========================================================================
    # Solve
    # ========================================================================

    println("Solving linear system...")

    K = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)

    t0 = time()
    sol_free = K \ reduced_rhs
    solve_time = time() - t0

    println("  Solve time: ", @sprintf("%.2f ms", solve_time * 1000))
    println()

    # ========================================================================
    # Reconstruct Fine-Scale Solution
    # ========================================================================

    println("Reconstructing fine-scale solution...")

    # Expand coarse solution to full DOFs (BCs = 0)
    coarse_solution = zeros(nnodes)
    for i = 1:nfree
        coarse_solution[free_to_global[i]] = sol_free[i]
    end

    t0 = time()
    fine_x, fine_y, fine_u = reconstruct_fine_solution(nx, ny, ratio, coarse_solution, phi_cpu, fem_mesh.cell_to_node)
    reconstruct_time = time() - t0

    println("  Fine mesh size: $(length(fine_u)) nodes")
    println("  Reconstruction time: ", @sprintf("%.2f ms", reconstruct_time * 1000))
    println()

    # ========================================================================
    # Output
    # ========================================================================

    println("Writing solutions to files...")

    # Use same file names as mfem_assembly.jl for easy comparison
    write_solution("mfem_solution_coarse.txt", fem_mesh.vertex_x, fem_mesh.vertex_y, coarse_solution)
    println("  Coarse solution written to: mfem_solution_coarse.txt")

    open("mfem_solution_fine.txt", "w") do io
        println(io, "# x y u")
        for i = 1:length(fine_u)
            @printf(io, "%.6f %.6f %.6e\n", fine_x[i], fine_y[i], fine_u[i])
        end
    end
    println("  Fine solution written to: mfem_solution_fine.txt")
    println()

    # ========================================================================
    # Performance summary (matching mfem_cuda.jl format)
    # ========================================================================

    # Group timings to match GPU version structure
    # "MFEM assembly" = precomp + CUDA kernel + GPU scatter
    # Do not include the transfer time to match the C++ code
    mfem_assembly_time = precomp_time + assembly_time + gpu_scatter_time

    println("="^70)
    println("Performance Summary")
    println("="^70)
    println("  Mesh generation:    ", @sprintf("%8.2f ms", mesh_time * 1000))
    println("  Connectivity/color: ", @sprintf("%8.2f ms", conn_time * 1000))
    println("  MFEM assembly:      ", @sprintf("%8.2f ms", mfem_assembly_time * 1000), "  (GPU)")
    println("    ├─ Pre-compute:   ", @sprintf("%8.2f ms (%.1f%%)", precomp_time * 1000, 100*precomp_time/mfem_assembly_time))
    println("    ├─ CUDA kernel:   ", @sprintf("%8.2f ms (%.1f%%)", assembly_time * 1000, 100*assembly_time/mfem_assembly_time))
    println("    ├─ GPU scatter:   ", @sprintf("%8.2f ms (%.1f%%)", gpu_scatter_time * 1000, 100*gpu_scatter_time/mfem_assembly_time))
    println("  GPU->CPU:           ", @sprintf("%8.2f ms", transfer_time * 1000))
    println("  Boundary conditions:", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Solve:              ", @sprintf("%8.2f ms", solve_time * 1000))
    println("  Reconstruction:     ", @sprintf("%8.2f ms", reconstruct_time * 1000))
    println("  " * "-"^68)
    total_workflow_time = mesh_time + conn_time + mfem_assembly_time + transfer_time + bc_time + solve_time + reconstruct_time
    println("  Total:              ", @sprintf("%8.2f ms", total_workflow_time * 1000))
    println()

    # Additional metrics for scaling studies
    println("="^70)
    println("Metrics for Scaling Studies")
    println("="^70)
    println("  Problem size:")
    println("    Coarse mesh:      $nx × $ny = ", nx*ny, " elements")
    println("    MFEM ratio:       $ratio")
    println("    Effective fine:   $(nx*ratio) × $(ny*ratio) = ", (nx*ratio)*(ny*ratio), " elements")
    println("    DOFs (coarse):    $nnodes")
    println("    Free DOFs:        $nfree")
    println()
    println("  Assembly breakdown (detailed):")
    println("    Pre-compute:      ", @sprintf("%8.2f ms (%.1f%%)", precomp_time * 1000, 100*precomp_time/mfem_assembly_time))
    println("    CUDA kernel:      ", @sprintf("%8.2f ms (%.1f%%)", assembly_time * 1000, 100*assembly_time/mfem_assembly_time))
    println("    GPU scatter:      ", @sprintf("%8.2f ms (%.1f%%)", gpu_scatter_time * 1000, 100*gpu_scatter_time/mfem_assembly_time))
    println("    GPU->CPU transfer:", @sprintf("%8.2f ms (%.1f%%)", transfer_time * 1000, 100*transfer_time/mfem_assembly_time))
    println("    Total assembly:   ", @sprintf("%8.2f ms", mfem_assembly_time * 1000))
    println()
    println("  Execution model:    GPU (CUDA kernels + coloring)")
    println("  Number of colors:   ", length(e2e_colors))
    println()
    println("  CUDA configuration:")
    println("    Blocks:           ", mesh.nel)
    println("    Threads/block:    256")
    println("    Shared memory:    $shmem_size bytes ($(shmem_size ÷ 1024) KB)")
    println("    Interior DOFs/element: ", precomp.nfree)
    println()
    println("  Throughput:")
    println("    Elements/sec:     ", @sprintf("%.0f", mesh.nel / assembly_time))
    println("    Fine elements/sec:", @sprintf("%.0f", mesh.nel * ratio * ratio / assembly_time))
    println()

    # Solution statistics
    println("Coarse solution statistics:")
    minval = minimum(sol_free)
    maxval = maximum(sol_free)
    avgval = sum(sol_free) / length(sol_free)
    println("  Min value: ", @sprintf("%.6e", minval))
    println("  Max value: ", @sprintf("%.6e", maxval))
    println("  Avg value: ", @sprintf("%.6e", avgval))
    println()

    println("Fine solution statistics:")
    minval_fine = minimum(fine_u)
    maxval_fine = maximum(fine_u)
    avgval_fine = sum(fine_u) / length(fine_u)
    println("  Min value: ", @sprintf("%.6e", minval_fine))
    println("  Max value: ", @sprintf("%.6e", maxval_fine))
    println("  Avg value: ", @sprintf("%.6e", avgval_fine))
    println()
end

# Run main program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
