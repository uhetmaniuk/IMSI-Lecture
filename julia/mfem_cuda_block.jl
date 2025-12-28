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
# Designed for ratio <= 16 to fit in shared memory (~35KB)
#

using CUDA
using LinearAlgebra
using SparseArrays
using Printf

# ============================================================================
# Constants and Configuration
# ============================================================================

# Maximum supported ratio (limited by shared memory)
const MAX_RATIO = 16
const MAX_FINE_NODES = (MAX_RATIO + 1)^2  # 289 for ratio=16
const MAX_NFREE = MAX_FINE_NODES - 4 * MAX_RATIO  # 225 for ratio=16
const MAX_NNZ_II = 9 * MAX_NFREE  # 2025 for ratio=16
const MAX_NBOUNDARY = 4 * MAX_RATIO  # 64 for ratio=16
const MAX_NNZ_B = 9 * MAX_NBOUNDARY  # 576 for ratio=16

# Number of coarse basis functions
const NUM_BASIS = 4
const NUM_BASIS_TO_SOLVE = 3

# PCG parameters
const PCG_TOL = 1e-10
const PCG_MAXITER = 500

# ============================================================================
# Gauss Quadrature (2x2 rule)
# ============================================================================

const GAUSS_PT = 1.0 / sqrt(3.0)
const GAUSS_PTS_XI = (-GAUSS_PT, GAUSS_PT, GAUSS_PT, -GAUSS_PT)
const GAUSS_PTS_ETA = (-GAUSS_PT, -GAUSS_PT, GAUSS_PT, GAUSS_PT)
const GAUSS_WTS = (1.0, 1.0, 1.0, 1.0)

# ============================================================================
# Device-side Q1 Shape Functions
# ============================================================================

"""
Evaluate Q1 shape functions at (xi, eta)
Returns tuple of 4 values
"""
@inline function shape_functions_device(xi::Float64, eta::Float64)
    N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
    N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
    N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
    N4 = 0.25 * (1.0 - xi) * (1.0 + eta)
    return (N1, N2, N3, N4)
end

"""
Evaluate Q1 shape function gradients at (xi, eta)
Returns tuple of 8 values: (dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN1_deta, dN2_deta, dN3_deta, dN4_deta)
"""
@inline function shape_gradients_device(xi::Float64, eta::Float64)
    # ∂N/∂ξ
    dN1_dxi = -0.25 * (1.0 - eta)
    dN2_dxi =  0.25 * (1.0 - eta)
    dN3_dxi =  0.25 * (1.0 + eta)
    dN4_dxi = -0.25 * (1.0 + eta)
    # ∂N/∂η
    dN1_deta = -0.25 * (1.0 - xi)
    dN2_deta = -0.25 * (1.0 + xi)
    dN3_deta =  0.25 * (1.0 + xi)
    dN4_deta =  0.25 * (1.0 - xi)
    return (dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN1_deta, dN2_deta, dN3_deta, dN4_deta)
end

# ============================================================================
# Block-Level Reduction Primitives
# ============================================================================

"""
Block-level sum reduction using shared memory
Returns the sum on all threads (broadcast result)
"""
@inline function block_reduce_sum(g, shmem_reduce, val::Float64, blocksize::Int32)
    tid = threadIdx().x

    # Store value in shared memory
    @inbounds shmem_reduce[tid] = val
    CG.sync(g)

    # Parallel tree reduction
    s = blocksize ÷ Int32(2)
    while s > Int32(0)
        if tid <= s
            @inbounds shmem_reduce[tid] += shmem_reduce[tid + s]
        end
        CG.sync(g)
        s ÷= Int32(2)
    end

    # Broadcast result to all threads
    @inbounds result = shmem_reduce[1]
    CG.sync(g)
    return result
end

"""
Block-level dot product: x'*y for vectors of length n
Each thread handles multiple elements via striding
"""
@inline function block_dot_product(g, shmem_reduce, x, y, n::Int32, blocksize::Int32)
    tid = threadIdx().x

    # Thread-local accumulation with striding
    local_sum = 0.0
    i = tid
    while i <= n
        @inbounds local_sum += x[i] * y[i]
        i += blocksize
    end

    return block_reduce_sum(g, shmem_reduce, local_sum, blocksize)
end

"""
Block-level sparse matrix-vector multiply: y = A*x
A is stored in CSC format (colptr, rowidx, nzval)
Each thread handles multiple columns via striding
"""
@inline function block_spmv!(g, y, n::Int32, colptr, rowidx, nzval, x, blocksize::Int32)
    tid = threadIdx().x

    # Each thread computes some rows (CSC: iterate over columns, accumulate to rows)
    # For symmetric matrix with CSC = CSR, we iterate over "rows"
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
# Block-Cooperative PCG Solver with Jacobi Preconditioner
# ============================================================================

"""
Block-cooperative PCG solver for a single RHS
All threads in block cooperate to solve A*x = b
Uses Jacobi (diagonal) preconditioning

Arguments:
- g: Thread block cooperative group
- x: Solution vector (in shared memory, modified in-place)
- n: System size
- colptr, rowidx, nzval: CSC sparse matrix arrays (in shared memory)
- b: Right-hand side (in shared memory)
- diag: Diagonal elements for Jacobi preconditioner (in shared memory)
- r, z, p, Ap: PCG workspace vectors (in shared memory)
- shmem_reduce: Shared memory for reductions
- blocksize: Number of threads in block
- tol: Convergence tolerance
- maxiter: Maximum iterations

Returns: Number of iterations (negative if not converged)
"""
@inline function block_pcg_jacobi!(g, x, n::Int32, colptr, rowidx, nzval, b, diag,
                                   r, z, p, Ap, shmem_reduce, blocksize::Int32,
                                   tol::Float64, maxiter::Int32)
    tid = threadIdx().x

    # Compute initial residual: r = b - A*x
    block_spmv!(g, r, n, colptr, rowidx, nzval, x, blocksize)
    i = tid
    while i <= n
        @inbounds r[i] = b[i] - r[i]
        i += blocksize
    end
    CG.sync(g)

    # Apply Jacobi preconditioner: z = D^{-1} * r
    # Initialize search direction: p = z
    i = tid
    while i <= n
        @inbounds d = diag[i]
        @inbounds z_val = (d != 0.0) ? r[i] / d : r[i]
        @inbounds z[i] = z_val
        @inbounds p[i] = z_val
        i += blocksize
    end
    CG.sync(g)

    # Compute initial rz = r'*z
    rz_old = block_dot_product(g, shmem_reduce, r, z, n, blocksize)

    # PCG iteration
    for iter = Int32(1):maxiter
        # Ap = A * p
        block_spmv!(g, Ap, n, colptr, rowidx, nzval, p, blocksize)

        # pAp = p' * Ap
        pAp = block_dot_product(g, shmem_reduce, p, Ap, n, blocksize)

        # Check for breakdown
        if pAp <= 0.0
            return -iter
        end

        # alpha = rz / pAp
        alpha = rz_old / pAp

        # x += alpha * p, r -= alpha * Ap
        i = tid
        while i <= n
            @inbounds x[i] += alpha * p[i]
            @inbounds r[i] -= alpha * Ap[i]
            i += blocksize
        end
        CG.sync(g)

        # Check convergence: ||r||^2
        rnorm_sq = block_dot_product(g, shmem_reduce, r, r, n, blocksize)

        if sqrt(rnorm_sq) < tol
            return iter
        end

        # z = D^{-1} * r (Jacobi preconditioner)
        i = tid
        while i <= n
            @inbounds d = diag[i]
            @inbounds z[i] = (d != 0.0) ? r[i] / d : r[i]
            i += blocksize
        end
        CG.sync(g)

        # rz_new = r' * z
        rz_new = block_dot_product(g, shmem_reduce, r, z, n, blocksize)

        # beta = rz_new / rz_old
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

    return -maxiter  # Did not converge
end

# ============================================================================
# Fine Element Assembly (Device Functions)
# ============================================================================

"""
Assemble a single fine element and return local stiffness matrix and RHS
Uses 2x2 Gauss quadrature
"""
@inline function assemble_fine_element_device(x1::Float64, y1::Float64, x2::Float64, y2::Float64,
                                               x3::Float64, y3::Float64, x4::Float64, y4::Float64,
                                               ax_val::Float64, ay_val::Float64, f_val::Float64)
    # Initialize local stiffness matrix and RHS (stored in registers)
    Ke11 = 0.0; Ke12 = 0.0; Ke13 = 0.0; Ke14 = 0.0
    Ke21 = 0.0; Ke22 = 0.0; Ke23 = 0.0; Ke24 = 0.0
    Ke31 = 0.0; Ke32 = 0.0; Ke33 = 0.0; Ke34 = 0.0
    Ke41 = 0.0; Ke42 = 0.0; Ke43 = 0.0; Ke44 = 0.0
    fe1 = 0.0; fe2 = 0.0; fe3 = 0.0; fe4 = 0.0

    # Loop over quadrature points
    for qp = 1:4
        xi = GAUSS_PTS_XI[qp]
        eta = GAUSS_PTS_ETA[qp]
        w = GAUSS_WTS[qp]

        # Shape functions
        N1, N2, N3, N4 = shape_functions_device(xi, eta)

        # Shape function gradients in reference coordinates
        dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN1_deta, dN2_deta, dN3_deta, dN4_deta = shape_gradients_device(xi, eta)

        # Jacobian matrix
        J11 = dN1_dxi * x1 + dN2_dxi * x2 + dN3_dxi * x3 + dN4_dxi * x4
        J12 = dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3 + dN4_deta * x4
        J21 = dN1_dxi * y1 + dN2_dxi * y2 + dN3_dxi * y3 + dN4_dxi * y4
        J22 = dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3 + dN4_deta * y4

        detJ = J11 * J22 - J12 * J21
        invdetJ = 1.0 / detJ

        # Inverse Jacobian
        invJ11 =  J22 * invdetJ
        invJ12 = -J12 * invdetJ
        invJ21 = -J21 * invdetJ
        invJ22 =  J11 * invdetJ

        # Gradients in physical coordinates
        dN1_dx = dN1_dxi * invJ11 + dN1_deta * invJ21
        dN1_dy = dN1_dxi * invJ12 + dN1_deta * invJ22
        dN2_dx = dN2_dxi * invJ11 + dN2_deta * invJ21
        dN2_dy = dN2_dxi * invJ12 + dN2_deta * invJ22
        dN3_dx = dN3_dxi * invJ11 + dN3_deta * invJ21
        dN3_dy = dN3_dxi * invJ12 + dN3_deta * invJ22
        dN4_dx = dN4_dxi * invJ11 + dN4_deta * invJ21
        dN4_dy = dN4_dxi * invJ12 + dN4_deta * invJ22

        # Integration weight
        detJ_w = detJ * w
        ax_detJ_w = ax_val * detJ_w
        ay_detJ_w = ay_val * detJ_w
        f_detJ_w = f_val * detJ_w

        # Accumulate stiffness matrix
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

        # Accumulate RHS
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

"""
Binary search to find position in CSC/CSR matrix
"""
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
# Main MFEM Element Kernel - One Block Per Coarse Element
# ============================================================================

"""
MFEM element assembly kernel - one block per coarse element

This kernel processes one coarse element per CUDA block. All threads in the
block cooperate to:
1. Build DOF mappings and sparsity patterns
2. Assemble fine elements in parallel
3. Solve for interior basis functions using block-cooperative PCG
4. Compute coarse element matrices via block reductions

Shared memory layout is dynamically computed based on ratio.
"""
function mfem_element_kernel!(Ke_out::CuDeviceMatrix{Float64},
                              fe_out::CuDeviceMatrix{Float64},
                              phi_out::CuDeviceArray{Float64, 3},
                              vertex_x::CuDeviceVector{Float64},
                              vertex_y::CuDeviceVector{Float64},
                              cell_to_node::CuDeviceMatrix{Int32},
                              ratio::Int32,
                              use_varying_coeffs::Bool)
    # Get cooperative group
    g = CG.this_thread_block()
    tid = threadIdx().x
    blocksize = Int32(blockDim().x)
    iel = Int32(blockIdx().x)  # One block per coarse element

    # Compute sizes
    numNodes = (ratio + Int32(1))^2
    nfree = numNodes - Int32(4) * ratio
    nboundary = Int32(4) * ratio
    max_nnz_ii = Int32(9) * nfree
    max_nnz_b = Int32(9) * nboundary

    # ========================================================================
    # Shared Memory Allocation
    # ========================================================================
    # Use dynamic shared memory with careful layout

    # Calculate offsets for Int32 section
    off_globalToFree = Int32(0)
    off_freeToGlobal = off_globalToFree + numNodes
    off_globalToBoundary = off_freeToGlobal + nfree
    off_boundaryToGlobal = off_globalToBoundary + numNodes
    off_colptr_ii = off_boundaryToGlobal + nboundary
    off_rowidx_ii = off_colptr_ii + nfree + Int32(1)
    off_matRowPtr_b = off_rowidx_ii + max_nnz_ii
    off_matColIdx_b = off_matRowPtr_b + nboundary + Int32(1)
    int32_total = off_matColIdx_b + max_nnz_b

    # Align to 8 bytes for Float64
    int32_total_aligned = ((int32_total + Int32(1)) ÷ Int32(2)) * Int32(2)

    # Float64 section offsets
    off_nzval_ii = Int32(0)
    off_matValues_b = off_nzval_ii + max_nnz_ii
    off_phi = off_matValues_b + max_nnz_b
    off_rhs_fine = off_phi + numNodes * Int32(4)
    off_btmp = off_rhs_fine + numNodes
    off_diag_ii = off_btmp + nfree * Int32(3)
    off_r = off_diag_ii + nfree
    off_z = off_r + nfree
    off_p = off_z + nfree
    off_Ap = off_p + nfree
    off_utmp = off_Ap + nfree
    off_Ke_coarse = off_utmp + nfree * Int32(3)
    off_fe_coarse = off_Ke_coarse + Int32(16)
    off_reduce = off_fe_coarse + Int32(4)

    # Allocate shared memory
    shmem_int32 = CuDynamicSharedArray(Int32, int32_total_aligned)
    shmem_float64 = CuDynamicSharedArray(Float64, off_reduce + blocksize, int32_total_aligned * sizeof(Int32))

    # Create views into shared memory
    globalToFree = @inbounds view(shmem_int32, (off_globalToFree + 1):(off_globalToFree + numNodes))
    freeToGlobal = @inbounds view(shmem_int32, (off_freeToGlobal + 1):(off_freeToGlobal + nfree))
    globalToBoundary = @inbounds view(shmem_int32, (off_globalToBoundary + 1):(off_globalToBoundary + numNodes))
    boundaryToGlobal = @inbounds view(shmem_int32, (off_boundaryToGlobal + 1):(off_boundaryToGlobal + nboundary))
    colptr_ii = @inbounds view(shmem_int32, (off_colptr_ii + 1):(off_colptr_ii + nfree + 1))
    rowidx_ii = @inbounds view(shmem_int32, (off_rowidx_ii + 1):(off_rowidx_ii + max_nnz_ii))
    matRowPtr_b = @inbounds view(shmem_int32, (off_matRowPtr_b + 1):(off_matRowPtr_b + nboundary + 1))
    matColIdx_b = @inbounds view(shmem_int32, (off_matColIdx_b + 1):(off_matColIdx_b + max_nnz_b))

    nzval_ii = @inbounds view(shmem_float64, (off_nzval_ii + 1):(off_nzval_ii + max_nnz_ii))
    matValues_b = @inbounds view(shmem_float64, (off_matValues_b + 1):(off_matValues_b + max_nnz_b))
    phi = @inbounds view(shmem_float64, (off_phi + 1):(off_phi + numNodes * 4))
    rhs_fine = @inbounds view(shmem_float64, (off_rhs_fine + 1):(off_rhs_fine + numNodes))
    btmp = @inbounds view(shmem_float64, (off_btmp + 1):(off_btmp + nfree * 3))
    diag_ii = @inbounds view(shmem_float64, (off_diag_ii + 1):(off_diag_ii + nfree))
    r = @inbounds view(shmem_float64, (off_r + 1):(off_r + nfree))
    z = @inbounds view(shmem_float64, (off_z + 1):(off_z + nfree))
    p = @inbounds view(shmem_float64, (off_p + 1):(off_p + nfree))
    Ap = @inbounds view(shmem_float64, (off_Ap + 1):(off_Ap + nfree))
    utmp = @inbounds view(shmem_float64, (off_utmp + 1):(off_utmp + nfree * 3))
    Ke_coarse = @inbounds view(shmem_float64, (off_Ke_coarse + 1):(off_Ke_coarse + 16))
    fe_coarse = @inbounds view(shmem_float64, (off_fe_coarse + 1):(off_fe_coarse + 4))
    shmem_reduce = @inbounds view(shmem_float64, (off_reduce + 1):(off_reduce + blocksize))

    # ========================================================================
    # Load Coarse Element Corner Coordinates
    # ========================================================================

    # Use first 8 positions of shmem_float64 temporarily for corner coords
    if tid <= Int32(4)
        @inbounds node = cell_to_node[iel, tid]
        @inbounds shmem_float64[2*tid - 1] = vertex_x[node]
        @inbounds shmem_float64[2*tid] = vertex_y[node]
    end
    CG.sync(g)

    @inbounds x_corner1 = shmem_float64[1]
    @inbounds y_corner1 = shmem_float64[2]
    @inbounds x_corner2 = shmem_float64[3]
    @inbounds y_corner2 = shmem_float64[4]
    @inbounds x_corner3 = shmem_float64[5]
    @inbounds y_corner3 = shmem_float64[6]
    @inbounds x_corner4 = shmem_float64[7]
    @inbounds y_corner4 = shmem_float64[8]

    # Compute fine element size
    hx = (x_corner2 - x_corner1) / Float64(ratio)
    hy = (y_corner4 - y_corner1) / Float64(ratio)

    # ========================================================================
    # Phase 1: Initialize Arrays and Build DOF Mappings
    # ========================================================================

    # Initialize arrays in parallel
    i = tid
    while i <= numNodes
        @inbounds globalToFree[i] = Int32(-1)
        @inbounds globalToBoundary[i] = Int32(-1)
        @inbounds phi[i] = 0.0
        @inbounds phi[i + numNodes] = 0.0
        @inbounds phi[i + 2*numNodes] = 0.0
        @inbounds phi[i + 3*numNodes] = 0.0
        @inbounds rhs_fine[i] = 0.0
        i += blocksize
    end
    i = tid
    while i <= max_nnz_ii
        @inbounds nzval_ii[i] = 0.0
        i += blocksize
    end
    i = tid
    while i <= max_nnz_b
        @inbounds matValues_b[i] = 0.0
        i += blocksize
    end
    i = tid
    while i <= nfree * Int32(3)
        @inbounds btmp[i] = 0.0
        i += blocksize
    end
    CG.sync(g)

    # Build DOF mappings (thread 1 only for consistency)
    if tid == Int32(1)
        nfree_count = Int32(0)
        nboundary_count = Int32(0)

        for iy = Int32(0):ratio
            for ix = Int32(0):ratio
                nodeID = ix + iy * (ratio + Int32(1)) + Int32(1)
                if ix > Int32(0) && ix < ratio && iy > Int32(0) && iy < ratio
                    nfree_count += Int32(1)
                    @inbounds freeToGlobal[nfree_count] = nodeID
                    @inbounds globalToFree[nodeID] = nfree_count
                else
                    nboundary_count += Int32(1)
                    @inbounds boundaryToGlobal[nboundary_count] = nodeID
                    @inbounds globalToBoundary[nodeID] = nboundary_count
                end
            end
        end

        # Build K_ii sparsity pattern
        @inbounds colptr_ii[1] = Int32(1)
        for i_free = Int32(1):nfree
            @inbounds iGlobal = freeToGlobal[i_free]
            ix = (iGlobal - Int32(1)) % (ratio + Int32(1))
            iy = (iGlobal - Int32(1)) ÷ (ratio + Int32(1))

            count = Int32(1)
            hasWest = (ix > Int32(1))
            hasEast = (ix < ratio - Int32(1))
            hasSouth = (iy > Int32(1))
            hasNorth = (iy < ratio - Int32(1))

            if hasSouth
                count += (hasWest ? Int32(1) : Int32(0)) + Int32(1) + (hasEast ? Int32(1) : Int32(0))
            end
            count += (hasWest ? Int32(1) : Int32(0)) + (hasEast ? Int32(1) : Int32(0))
            if hasNorth
                count += (hasWest ? Int32(1) : Int32(0)) + Int32(1) + (hasEast ? Int32(1) : Int32(0))
            end

            @inbounds colptr_ii[i_free + 1] = colptr_ii[i_free] + count
        end

        # Fill K_ii column indices
        for i_free = Int32(1):nfree
            @inbounds iGlobal = freeToGlobal[i_free]
            ix = (iGlobal - Int32(1)) % (ratio + Int32(1))
            iy = (iGlobal - Int32(1)) ÷ (ratio + Int32(1))
            @inbounds offset = colptr_ii[i_free]

            if iy > Int32(1)
                if ix > Int32(1)
                    jGlobal = iGlobal - Int32(1) - (ratio + Int32(1))
                    @inbounds jFree = globalToFree[jGlobal]
                    if jFree != Int32(-1)
                        @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                    end
                end
                jGlobal = iGlobal - (ratio + Int32(1))
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != Int32(-1)
                    @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                end
                if ix < ratio - Int32(1)
                    jGlobal = iGlobal + Int32(1) - (ratio + Int32(1))
                    @inbounds jFree = globalToFree[jGlobal]
                    if jFree != Int32(-1)
                        @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                    end
                end
            end
            if ix > Int32(1)
                jGlobal = iGlobal - Int32(1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != Int32(-1)
                    @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                end
            end
            @inbounds rowidx_ii[offset] = i_free; offset += Int32(1)
            if ix < ratio - Int32(1)
                jGlobal = iGlobal + Int32(1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != Int32(-1)
                    @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                end
            end
            if iy < ratio - Int32(1)
                if ix > Int32(1)
                    jGlobal = iGlobal - Int32(1) + (ratio + Int32(1))
                    @inbounds jFree = globalToFree[jGlobal]
                    if jFree != Int32(-1)
                        @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                    end
                end
                jGlobal = iGlobal + (ratio + Int32(1))
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != Int32(-1)
                    @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                end
                if ix < ratio - Int32(1)
                    jGlobal = iGlobal + Int32(1) + (ratio + Int32(1))
                    @inbounds jFree = globalToFree[jGlobal]
                    if jFree != Int32(-1)
                        @inbounds rowidx_ii[offset] = jFree; offset += Int32(1)
                    end
                end
            end
        end

        # Build K_b sparsity pattern
        @inbounds matRowPtr_b[1] = Int32(1)
        for i_bnd = Int32(1):nboundary
            @inbounds iGlobal = boundaryToGlobal[i_bnd]
            ix = (iGlobal - Int32(1)) % (ratio + Int32(1))
            iy = (iGlobal - Int32(1)) ÷ (ratio + Int32(1))

            count = Int32(1)
            hasWest = (ix > Int32(0))
            hasEast = (ix < ratio)

            count += (hasWest ? Int32(1) : Int32(0)) + (hasEast ? Int32(1) : Int32(0))
            if iy > Int32(0)
                count += Int32(1) + (hasWest ? Int32(1) : Int32(0)) + (hasEast ? Int32(1) : Int32(0))
            end
            if iy < ratio
                count += Int32(1) + (hasWest ? Int32(1) : Int32(0)) + (hasEast ? Int32(1) : Int32(0))
            end

            @inbounds matRowPtr_b[i_bnd + 1] = matRowPtr_b[i_bnd] + count
        end

        # Fill K_b column indices
        for i_bnd = Int32(1):nboundary
            @inbounds iGlobal = boundaryToGlobal[i_bnd]
            ix = (iGlobal - Int32(1)) % (ratio + Int32(1))
            iy = (iGlobal - Int32(1)) ÷ (ratio + Int32(1))
            @inbounds offset = matRowPtr_b[i_bnd]

            if iy > Int32(0)
                if ix > Int32(0)
                    @inbounds matColIdx_b[offset] = iGlobal - Int32(1) - (ratio + Int32(1)); offset += Int32(1)
                end
                @inbounds matColIdx_b[offset] = iGlobal - (ratio + Int32(1)); offset += Int32(1)
                if ix < ratio
                    @inbounds matColIdx_b[offset] = iGlobal + Int32(1) - (ratio + Int32(1)); offset += Int32(1)
                end
            end
            if ix > Int32(0)
                @inbounds matColIdx_b[offset] = iGlobal - Int32(1); offset += Int32(1)
            end
            @inbounds matColIdx_b[offset] = iGlobal; offset += Int32(1)
            if ix < ratio
                @inbounds matColIdx_b[offset] = iGlobal + Int32(1); offset += Int32(1)
            end
            if iy < ratio
                if ix > Int32(0)
                    @inbounds matColIdx_b[offset] = iGlobal - Int32(1) + (ratio + Int32(1)); offset += Int32(1)
                end
                @inbounds matColIdx_b[offset] = iGlobal + (ratio + Int32(1)); offset += Int32(1)
                if ix < ratio
                    @inbounds matColIdx_b[offset] = iGlobal + Int32(1) + (ratio + Int32(1)); offset += Int32(1)
                end
            end
        end

        # Initialize boundary phi values
        @inbounds phi[1] = 1.0
        @inbounds phi[ratio + 1 + numNodes] = 1.0
        @inbounds phi[numNodes + 2*numNodes] = 1.0
        @inbounds phi[(ratio + 1) * ratio + 1 + 3*numNodes] = 1.0

        for is = Int32(1):(ratio - Int32(1))
            s = Float64(is) / Float64(ratio)
            nodeID = is * (ratio + Int32(1)) + Int32(1)
            @inbounds phi[nodeID] = 1.0 - s
            @inbounds phi[nodeID + 3*numNodes] = s
            nodeID = ratio + is * (ratio + Int32(1)) + Int32(1)
            @inbounds phi[nodeID + numNodes] = 1.0 - s
            @inbounds phi[nodeID + 2*numNodes] = s
            nodeID = is + Int32(1)
            @inbounds phi[nodeID] = 1.0 - s
            @inbounds phi[nodeID + numNodes] = s
            nodeID = is + ratio * (ratio + Int32(1)) + Int32(1)
            @inbounds phi[nodeID + 3*numNodes] = 1.0 - s
            @inbounds phi[nodeID + 2*numNodes] = s
        end
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

        xc = 0.25 * (x1 + x2 + x3 + x4)
        yc = 0.25 * (y1 + y2 + y3 + y4)

        ax_val = use_varying_coeffs ? (1.0 + 100.0 * cos(150.0 * xc)^2 * sin(150.0 * yc)^2) : 1.0
        ay_val = ax_val
        f_val = use_varying_coeffs ? sin(xc) * sin(yc) : 2.0 * 3.141592653589793^2 * sin(3.141592653589793 * xc) * sin(3.141592653589793 * yc)

        Ke11, Ke12, Ke13, Ke14, Ke21, Ke22, Ke23, Ke24, Ke31, Ke32, Ke33, Ke34, Ke41, Ke42, Ke43, Ke44, fe1, fe2, fe3, fe4 =
            assemble_fine_element_device(x1, y1, x2, y2, x3, y3, x4, y4, ax_val, ay_val, f_val)

        node1 = ix + iy * (ratio + Int32(1)) + Int32(1)
        node2 = node1 + Int32(1)
        node3 = node2 + (ratio + Int32(1))
        node4 = node1 + (ratio + Int32(1))

        Ke = (Ke11, Ke12, Ke13, Ke14, Ke21, Ke22, Ke23, Ke24, Ke31, Ke32, Ke33, Ke34, Ke41, Ke42, Ke43, Ke44)
        fe = (fe1, fe2, fe3, fe4)
        nodes = (node1, node2, node3, node4)

        for i = 1:4
            @inbounds gi = nodes[i]
            @inbounds CUDA.@atomic rhs_fine[gi] += fe[i]
        end

        for i = 1:4
            @inbounds iGlobal = nodes[i]
            @inbounds iFree = globalToFree[iGlobal]

            if iFree != Int32(-1)
                for j = 1:4
                    @inbounds jGlobal = nodes[j]
                    @inbounds jFree = globalToFree[jGlobal]
                    @inbounds k_val = Ke[(i-1)*4 + j]

                    if jFree != Int32(-1)
                        k = find_matrix_position_device(colptr_ii, rowidx_ii, iFree, jFree)
                        if k != Int32(-1)
                            @inbounds CUDA.@atomic nzval_ii[k] += k_val
                        end
                    else
                        for ir = 1:3
                            @inbounds phi_val = phi[jGlobal + (ir-1)*numNodes]
                            @inbounds CUDA.@atomic btmp[iFree + (ir-1)*nfree] -= k_val * phi_val
                        end
                    end
                end
            else
                @inbounds iBoundary = globalToBoundary[iGlobal]
                for j = 1:4
                    @inbounds jGlobal = nodes[j]
                    @inbounds k_val = Ke[(i-1)*4 + j]
                    k = find_matrix_position_device(matRowPtr_b, matColIdx_b, iBoundary, jGlobal)
                    if k != Int32(-1)
                        @inbounds CUDA.@atomic matValues_b[k] += k_val
                    end
                end
            end
        end

        ielem += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 3: Extract Diagonal and Build Initial Guess
    # ========================================================================

    i = tid
    while i <= nfree
        @inbounds diag_ii[i] = 0.0
        @inbounds for k = colptr_ii[i]:(colptr_ii[i + 1] - Int32(1))
            if rowidx_ii[k] == i
                diag_ii[i] = nzval_ii[k]
                break
            end
        end
        i += blocksize
    end
    CG.sync(g)

    i = tid
    while i <= nfree
        @inbounds gi = freeToGlobal[i]
        iy_fine = (gi - Int32(1)) ÷ (ratio + Int32(1))
        ix_fine = (gi - Int32(1)) % (ratio + Int32(1))
        xi = Float64(ix_fine) / Float64(ratio)
        eta = Float64(iy_fine) / Float64(ratio)

        N1 = (1.0 - xi) * (1.0 - eta)
        N2 = xi * (1.0 - eta)
        N3 = xi * eta

        @inbounds utmp[i] = N1
        @inbounds utmp[i + nfree] = N2
        @inbounds utmp[i + 2*nfree] = N3
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 4: Solve for Interior Basis Functions with Block-Cooperative PCG
    # ========================================================================

    for ir = Int32(1):Int32(3)
        x_offset = (ir - Int32(1)) * nfree
        b_offset = (ir - Int32(1)) * nfree

        x_view = @inbounds view(utmp, (x_offset + 1):(x_offset + nfree))
        b_view = @inbounds view(btmp, (b_offset + 1):(b_offset + nfree))

        block_pcg_jacobi!(g, x_view, nfree, colptr_ii, rowidx_ii, nzval_ii, b_view, diag_ii,
                          r, z, p, Ap, shmem_reduce, blocksize,
                          PCG_TOL, Int32(PCG_MAXITER))
    end
    CG.sync(g)

    # ========================================================================
    # Phase 5: Reconstruct Full Basis Functions
    # ========================================================================

    i = tid
    while i <= nfree
        @inbounds gi = freeToGlobal[i]
        sum_val = 0.0
        for ir = 1:3
            @inbounds val = utmp[i + (ir-1)*nfree]
            @inbounds phi[gi + (ir-1)*numNodes] = val
            sum_val += val
        end
        @inbounds phi[gi + 3*numNodes] = 1.0 - sum_val
        i += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 6: Compute Coarse Element Matrices
    # ========================================================================

    if tid <= Int32(16)
        @inbounds Ke_coarse[tid] = 0.0
    end
    if tid <= Int32(4)
        @inbounds fe_coarse[tid] = 0.0
    end
    CG.sync(g)

    for ir = Int32(1):Int32(4)
        dot_val = block_dot_product(g, shmem_reduce,
                                    view(phi, ((ir-1)*numNodes + 1):(ir*numNodes)),
                                    rhs_fine, numNodes, blocksize)
        if tid == Int32(1)
            @inbounds fe_coarse[ir] = dot_val
        end
        CG.sync(g)
    end

    i_bnd = tid
    while i_bnd <= nboundary
        @inbounds iGlobal = boundaryToGlobal[i_bnd]
        @inbounds rowBegin = matRowPtr_b[i_bnd]
        @inbounds rowEnd = matRowPtr_b[i_bnd + 1] - Int32(1)

        for k = rowBegin:rowEnd
            @inbounds jGlobal = matColIdx_b[k]
            @inbounds k_val = matValues_b[k]

            for ir = 1:4
                @inbounds phi_i = phi[iGlobal + (ir-1)*numNodes]
                phi_i_k = phi_i * k_val
                for jr = 1:4
                    @inbounds phi_j = phi[jGlobal + (jr-1)*numNodes]
                    @inbounds CUDA.@atomic Ke_coarse[(ir-1)*4 + jr] += phi_i_k * phi_j
                end
            end
        end
        i_bnd += blocksize
    end
    CG.sync(g)

    # ========================================================================
    # Phase 7: Write Output
    # ========================================================================

    if tid <= Int32(16)
        @inbounds Ke_out[iel, tid] = Ke_coarse[tid]
    end
    if tid <= Int32(4)
        @inbounds fe_out[iel, tid] = fe_coarse[tid]
    end

    i = tid
    while i <= numNodes
        for ir = 1:4
            @inbounds phi_out[iel, i, ir] = phi[i + (ir-1)*numNodes]
        end
        i += blocksize
    end

    return nothing
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

function calculate_shmem_size(ratio::Int, blocksize::Int)
    numNodes = (ratio + 1)^2
    nfree = numNodes - 4 * ratio
    nboundary = 4 * ratio
    max_nnz_ii = 9 * nfree
    max_nnz_b = 9 * nboundary

    int32_count = (2 * numNodes + nfree + nboundary +
                   nfree + 1 + max_nnz_ii +
                   nboundary + 1 + max_nnz_b)
    int32_aligned = ((int32_count + 1) ÷ 2) * 2

    float64_count = (max_nnz_ii + max_nnz_b +
                     numNodes * 4 + numNodes +
                     nfree * 3 + nfree +
                     4 * nfree + nfree * 3 +
                     16 + 4 + blocksize)

    return int32_aligned * sizeof(Int32) + float64_count * sizeof(Float64)
end

function assemble_mfem_cuda_block(mesh::CudaMesh, ratio::Int;
                                   blocksize::Int=256,
                                   use_varying_coeffs::Bool=true)
    nel = mesh.nel
    numNodes = (ratio + 1)^2

    if ratio > MAX_RATIO
        error("Ratio $ratio exceeds maximum supported ratio $MAX_RATIO")
    end

    Ke_out = CUDA.zeros(Float64, nel, 16)
    fe_out = CUDA.zeros(Float64, nel, 4)
    phi_out = CUDA.zeros(Float64, nel, numNodes, 4)

    shmem_size = calculate_shmem_size(ratio, blocksize)

    dev = CUDA.device()
    max_shmem = CUDA.attribute(dev, CUDA.CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN)
    if shmem_size > max_shmem
        error("Required shared memory ($shmem_size bytes) exceeds device limit ($max_shmem bytes)")
    end

    @cuda threads=blocksize blocks=nel shmem=shmem_size mfem_element_kernel!(
        Ke_out, fe_out, phi_out,
        mesh.vertex_x, mesh.vertex_y, mesh.cell_to_node,
        Int32(ratio), use_varying_coeffs
    )

    return Ke_out, fe_out, phi_out
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
    println("CUDA MFEM Assembly Test (One Block Per Coarse Element)")
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

    vertex_x, vertex_y, cell_to_node = generate_test_mesh(nx, ny)
    mesh = CudaMesh(vertex_x, vertex_y, cell_to_node)

    shmem_size = calculate_shmem_size(ratio, 256)
    println("Shared memory per block: $(shmem_size ÷ 1024) KB")
    println()

    println("Running warm-up...")
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, ratio)
    CUDA.synchronize()

    println("Running timed assembly...")
    t_start = time()
    Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, ratio)
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
    println("CUDA MFEM Assembly Benchmark (One Block Per Element)")
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

    vertex_x, vertex_y, cell_to_node = generate_test_mesh(nx, ny)
    mesh = CudaMesh(vertex_x, vertex_y, cell_to_node)

    for _ = 1:3
        Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, ratio)
        CUDA.synchronize()
    end

    times = Float64[]
    for _ = 1:nruns
        t = CUDA.@elapsed begin
            Ke_out, fe_out, phi_out = assemble_mfem_cuda_block(mesh, ratio)
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

if abspath(PROGRAM_FILE) == @__FILE__
    test_mfem_cuda_block()
    println()
    benchmark_mfem_cuda_block(16, 16, 8)
end
