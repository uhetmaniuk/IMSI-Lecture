#!/usr/bin/env julia
#
# @file mfem_assembly.jl
# @brief Julia implementation of 2D MFEM (Multiscale Finite Element Method)
#
# Solves: -∇·(α∇u) = f on a 2D rectangular domain
# with homogeneous Dirichlet boundary conditions
#
# Features:
# - 2D MFEM elements with local fine-grid resolution
# - Varying coefficients ax, ay, f
# - Thread-parallel assembly with graph coloring
# - Built-in sparse solver (\)
# - Static condensation for coarse-scale system
#

using CUDA
using CUDA.CUSPARSE
using Krylov
import Krylov: CgWorkspace, cg!
using SparseArrays
using LinearAlgebra
using Printf

# Fix for thread scaling: disable nested threading in BLAS
# This prevents oversubscription when using Threads.@threads
BLAS.set_num_threads(1)

# Load shared FEM base kernels
include("FEMBase.jl")
using .FEMBase

# Load PCG solver
include("PCG.jl")
using .PCG

# Load MFEM workspace
include("MFEMWorkspace.jl")
using .MFEMWorkspace


# ========================================================================
# CUDA Kernels
# ========================================================================

# Device functions for material coefficients and forcing term
# Set USE_VARYING_COEFFICIENTS to true to test with non-constant diffusion
const USE_VARYING_COEFFICIENTS = false

@inline function ax_func(x, y)
    if USE_VARYING_COEFFICIENTS
        # Varying coefficient: 1 + 0.5*sin(2πx)*sin(2πy)
        return 1.0 + 0.5 * sin(2π * x) * sin(2π * y)
    else
        return 1.0  # Constant diffusion coefficient
    end
end

@inline function ay_func(x, y)
    if USE_VARYING_COEFFICIENTS
        # Varying coefficient: 1 + 0.5*cos(2πx)*cos(2πy)
        return 1.0 + 0.5 * cos(2π * x) * cos(2π * y)
    else
        return 1.0  # Constant diffusion coefficient
    end
end

@inline function f_func(x, y)
    # Manufactured solution: u = sin(π*x) * sin(π*y)
    # Note: For varying coefficients, this RHS is no longer exact
    return 2.0 * π^2 * sin(π * x) * sin(π * y)
end

# Helper function to get Ke value with constant tuple indices (GPU-safe)
@inline function get_Ke_val(Ke, i, j)
    # Use if-else chain to avoid dynamic tuple indexing
    idx = (i - 1) * 4 + j
    if idx == 1; return Ke[1]
    elseif idx == 2; return Ke[2]
    elseif idx == 3; return Ke[3]
    elseif idx == 4; return Ke[4]
    elseif idx == 5; return Ke[5]
    elseif idx == 6; return Ke[6]
    elseif idx == 7; return Ke[7]
    elseif idx == 8; return Ke[8]
    elseif idx == 9; return Ke[9]
    elseif idx == 10; return Ke[10]
    elseif idx == 11; return Ke[11]
    elseif idx == 12; return Ke[12]
    elseif idx == 13; return Ke[13]
    elseif idx == 14; return Ke[14]
    elseif idx == 15; return Ke[15]
    else return Ke[16]
    end
end

# Helper function to get node from nodeList tuple (GPU-safe)
@inline function get_node(nodes, i)
    if i == 1; return nodes[1]
    elseif i == 2; return nodes[2]
    elseif i == 3; return nodes[3]
    else return nodes[4]
    end
end

# CUDA kernel to extract diagonal elements from CSC sparse matrix
function extract_diagonal_kernel!(d_diag, d_colptr, d_rowidx, d_nzval, n)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= n
        # Search for diagonal entry in column i
        @inbounds col_start = d_colptr[i]
        @inbounds col_end = d_colptr[i + 1] - 1
        for k = col_start:col_end
            @inbounds if d_rowidx[k] == i
                @inbounds d_diag[i] = d_nzval[k]
                break
            end
        end
    end
    return nothing
end

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

# CUDA kernel to assemble K_ii, btmp (RHS), and rhs_fine for all
# coarse elements
function assemble_Kii_kernel!(
    d_valK_ii, d_btmp, d_rhs_fine, d_phi, d_colptr, d_rowidx,
    d_globalToFree, d_valK_b, d_rowptr_b, d_colidx_b, d_globalToBoundary,
    nx, ny, ratio, log2_ratio, nfree, nnz_ii, nboundary, nnz_b, numVectorsToSolve,
    numNodes
)
    # One thread per coarse element
    iel = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if iel <= nx * ny
        # Coarse element indices (0-based)
        iel_0 = iel - 1
        ix_coarse = iel_0 % nx
        iy_coarse = iel_0 ÷ nx

        # Coarse element corners in [0,1]^2 domain
        hx_coarse = 1.0 / nx
        hy_coarse = 1.0 / ny
        x_c0 = Float64(ix_coarse) * hx_coarse
        y_c0 = Float64(iy_coarse) * hy_coarse

        # Fine element size within this coarse element
        # Use ldexp for exact division by power of 2
        hx = ldexp(hx_coarse, -log2_ratio)
        hy = ldexp(hy_coarse, -log2_ratio)

        # Offset for this block's K_ii and K_b in the global arrays
        offset_ii = (iel - 1) * nnz_ii
        offset_b = (iel - 1) * nnz_b
        offset_nodes = (iel - 1) * numNodes

        # Loop over fine elements within this coarse element
        for iy_fine = 0:ratio-1
            for ix_fine = 0:ratio-1
                # Fine element corners
                x1 = x_c0 + ix_fine * hx
                y1 = y_c0 + iy_fine * hy

                # Assemble Q1 element stiffness matrix (2x2 Gauss quad)
                # For axis-aligned rectangular elements, Jacobian is const
                J11 = 0.5 * hx
                J22 = 0.5 * hy
                detJ = J11 * J22
                invJ11 = 1.0 / J11
                invJ22 = 1.0 / J22

                # Gauss quadrature: 2x2 rule
                gp = 1.0 / sqrt(3.0)

                # Initialize element stiffness and RHS (explicit scalars)
                Ke11 = 0.0; Ke12 = 0.0; Ke13 = 0.0; Ke14 = 0.0
                Ke21 = 0.0; Ke22 = 0.0; Ke23 = 0.0; Ke24 = 0.0
                Ke31 = 0.0; Ke32 = 0.0; Ke33 = 0.0; Ke34 = 0.0
                Ke41 = 0.0; Ke42 = 0.0; Ke43 = 0.0; Ke44 = 0.0
                fe1 = 0.0; fe2 = 0.0; fe3 = 0.0; fe4 = 0.0

                # Loop over 4 Gauss points (unrolled)
                for qp = 1:4
                    xi = qp == 1 ? -gp : (qp == 2 ? gp : (qp == 3 ? gp : -gp))
                    eta = qp <= 2 ? -gp : gp

                    # Q1 shape functions at Gauss point
                    N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
                    N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
                    N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
                    N4 = 0.25 * (1.0 - xi) * (1.0 + eta)

                    # Physical coordinates at Gauss point
                    xq = x1 * N1 + (x1 + hx) * N2 + (x1 + hx) * N3 + x1 * N4
                    yq = y1 * N1 + y1 * N2 + (y1 + hy) * N3 + (y1 + hy) * N4

                    # Evaluate coefficients at Gauss point
                    ax_val = ax_func(xq, yq)
                    ay_val = ay_func(xq, yq)
                    f_val = f_func(xq, yq)

                    # Q1 shape function gradients in reference coords
                    dN1_dxi = -0.25 * (1.0 - eta)
                    dN2_dxi =  0.25 * (1.0 - eta)
                    dN3_dxi =  0.25 * (1.0 + eta)
                    dN4_dxi = -0.25 * (1.0 + eta)

                    dN1_deta = -0.25 * (1.0 - xi)
                    dN2_deta = -0.25 * (1.0 + xi)
                    dN3_deta =  0.25 * (1.0 + xi)
                    dN4_deta =  0.25 * (1.0 - xi)

                    # Physical gradients
                    dN1_dx = dN1_dxi * invJ11; dN1_dy = dN1_deta * invJ22
                    dN2_dx = dN2_dxi * invJ11; dN2_dy = dN2_deta * invJ22
                    dN3_dx = dN3_dxi * invJ11; dN3_dy = dN3_deta * invJ22
                    dN4_dx = dN4_dxi * invJ11; dN4_dy = dN4_deta * invJ22

                    # Quadrature weight * detJ (weight = 1 for 2x2 Gauss)
                    w_detJ = detJ

                    # Accumulate stiffness (all 16 entries)
                    Ke11 += (ax_val * dN1_dx * dN1_dx + ay_val * dN1_dy * dN1_dy) * w_detJ
                    Ke12 += (ax_val * dN1_dx * dN2_dx + ay_val * dN1_dy * dN2_dy) * w_detJ
                    Ke13 += (ax_val * dN1_dx * dN3_dx + ay_val * dN1_dy * dN3_dy) * w_detJ
                    Ke14 += (ax_val * dN1_dx * dN4_dx + ay_val * dN1_dy * dN4_dy) * w_detJ

                    Ke21 += (ax_val * dN2_dx * dN1_dx + ay_val * dN2_dy * dN1_dy) * w_detJ
                    Ke22 += (ax_val * dN2_dx * dN2_dx + ay_val * dN2_dy * dN2_dy) * w_detJ
                    Ke23 += (ax_val * dN2_dx * dN3_dx + ay_val * dN2_dy * dN3_dy) * w_detJ
                    Ke24 += (ax_val * dN2_dx * dN4_dx + ay_val * dN2_dy * dN4_dy) * w_detJ

                    Ke31 += (ax_val * dN3_dx * dN1_dx + ay_val * dN3_dy * dN1_dy) * w_detJ
                    Ke32 += (ax_val * dN3_dx * dN2_dx + ay_val * dN3_dy * dN2_dy) * w_detJ
                    Ke33 += (ax_val * dN3_dx * dN3_dx + ay_val * dN3_dy * dN3_dy) * w_detJ
                    Ke34 += (ax_val * dN3_dx * dN4_dx + ay_val * dN3_dy * dN4_dy) * w_detJ

                    Ke41 += (ax_val * dN4_dx * dN1_dx + ay_val * dN4_dy * dN1_dy) * w_detJ
                    Ke42 += (ax_val * dN4_dx * dN2_dx + ay_val * dN4_dy * dN2_dy) * w_detJ
                    Ke43 += (ax_val * dN4_dx * dN3_dx + ay_val * dN4_dy * dN3_dy) * w_detJ
                    Ke44 += (ax_val * dN4_dx * dN4_dx + ay_val * dN4_dy * dN4_dy) * w_detJ

                    # Accumulate RHS
                    fe1 += N1 * f_val * w_detJ
                    fe2 += N2 * f_val * w_detJ
                    fe3 += N3 * f_val * w_detJ
                    fe4 += N4 * f_val * w_detJ
                end

                # Store Ke in row-major tuple for later access
                Ke_fine = (Ke11, Ke12, Ke13, Ke14,
                           Ke21, Ke22, Ke23, Ke24,
                           Ke31, Ke32, Ke33, Ke34,
                           Ke41, Ke42, Ke43, Ke44)

                # Fine element node list (local numbering within coarse elem)
                node1 = ix_fine + iy_fine * (ratio + 1) + 1
                node2 = ix_fine + 1 + iy_fine * (ratio + 1) + 1
                node3 = ix_fine + 1 + (iy_fine + 1) * (ratio + 1) + 1
                node4 = ix_fine + (iy_fine + 1) * (ratio + 1) + 1

                # Scatter RHS to global d_rhs_fine (explicit to avoid dynamic dispatch)
                offset_rhs = (iel - 1) * numNodes
                CUDA.@atomic d_rhs_fine[offset_rhs + node1] += fe1
                CUDA.@atomic d_rhs_fine[offset_rhs + node2] += fe2
                CUDA.@atomic d_rhs_fine[offset_rhs + node3] += fe3
                CUDA.@atomic d_rhs_fine[offset_rhs + node4] += fe4

                # Create nodeList tuple for scatter loops
                nodeList = (node1, node2, node3, node4)

                # Scatter to K_ii and K_b
                # Block offsets for global indexing into block-diagonal matrices
                block_offset_free = (iel - 1) * nfree
                block_offset_boundary = (iel - 1) * nboundary

                for i = 1:4
                    iGlobal = get_node(nodeList, i)
                    iFree_local = d_globalToFree[iGlobal]
                    iBoundary_local = d_globalToBoundary[iGlobal]

                    # Scatter to K_ii (only for interior nodes)
                    if iFree_local != -1
                        # Global index into block-diagonal K_ii
                        iFree_global = block_offset_free + iFree_local

                        for j = 1:4
                            jGlobal = get_node(nodeList, j)
                            jFree_local = d_globalToFree[jGlobal]
                            k_val = get_Ke_val(Ke_fine, i, j)

                            if jFree_local != -1
                                # Interior-interior: add to K_ii
                                # Use global column index for colptr lookup
                                col_start = d_colptr[iFree_global]
                                col_end = d_colptr[iFree_global + 1] - 1

                                # Global row index to search for
                                jFree_global = block_offset_free + jFree_local

                                # Binary search for global jFree in rowidx
                                k_pos = binary_search_sparse(
                                    d_rowidx, jFree_global, col_start, col_end
                                )
                                if k_pos != -1
                                    # k_pos is already global, no offset needed
                                    CUDA.@atomic d_valK_ii[k_pos] += k_val
                                end
                            else
                                # Interior-boundary: add to btmp
                                # (RHS for static condensation)
                                # Offset for this block's btmp
                                offset_btmp = (iel - 1) * nfree

                                # Unroll loop over 3 basis functions
                                phi_val1 = d_phi[offset_nodes + jGlobal, 1]
                                phi_val2 = d_phi[offset_nodes + jGlobal, 2]
                                phi_val3 = d_phi[offset_nodes + jGlobal, 3]
                                btmp_idx = offset_btmp + iFree_local
                                CUDA.@atomic d_btmp[btmp_idx, 1] -= k_val * phi_val1
                                CUDA.@atomic d_btmp[btmp_idx, 2] -= k_val * phi_val2
                                CUDA.@atomic d_btmp[btmp_idx, 3] -= k_val * phi_val3
                            end
                        end
                    end

                    # Scatter to K_b (only for boundary nodes)
                    if iBoundary_local != -1
                        # Global row index into block-diagonal K_b
                        iBoundary_global = block_offset_boundary + iBoundary_local

                        for j = 1:4
                            jGlobal = get_node(nodeList, j)
                            jGlobal_offset = offset_nodes + jGlobal
                            k_val = get_Ke_val(Ke_fine, i, j)

                            # Boundary-all: add to K_b
                            # Use global row index for rowptr lookup
                            row_start = d_rowptr_b[iBoundary_global]
                            row_end = d_rowptr_b[iBoundary_global + 1] - 1

                            # Binary search for jGlobal_offset in colidx
                            k_pos = binary_search_sparse(
                                d_colidx_b, jGlobal_offset,
                                row_start, row_end
                            )
                            if k_pos != -1
                                # k_pos is already global, no offset needed
                                CUDA.@atomic d_valK_b[k_pos] += k_val
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# CUDA kernel to build block diagonal sparsity pattern
# (optimized for full thread utilization)
function build_block_diagonal_sparsity_kernel!(
    d_colptr, d_rowidx, d_colptrLocal, d_rowidxLocal, nb, nfree, nnz_ii
)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    # Each thread handles one colptr entry (1-based Julia indexing)
    if idx < nb * (nfree + 1)
        k = idx ÷ (nfree + 1)  # Block index
        j = idx % (nfree + 1)   # Local column index (0 to nfree)
        offset_col = k * nfree
        offset_nnz = k * nnz_ii
        # +1 for Julia 1-based indexing
        @inbounds d_colptr[offset_col + j + 1] = (
            d_colptrLocal[j + 1] + offset_nnz
        )
    end

    # Each thread handles one rowidx entry (offset to avoid overlap)
    idx_row = idx - nb * (nfree + 1)
    if idx_row >= 0 && idx_row < nb * nnz_ii
        k = idx_row ÷ nnz_ii
        i = idx_row % nnz_ii
        offset_nnz = k * nnz_ii
        @inbounds d_rowidx[offset_nnz + i + 1] = (
            d_rowidxLocal[i + 1] + k * nfree
        )
    end

    return nothing
end

# CUDA kernel to build block diagonal CSR sparsity pattern for K_b
# (optimized for full thread utilization)
function build_block_diagonal_csr_sparsity_kernel!(
    d_rowptr, d_colidx, d_rowptrLocal, d_colidxLocal,
    nb, nboundary, nnz_b, numNodes
)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    # Each thread handles one rowptr entry (1-based Julia indexing)
    if idx < nb * (nboundary + 1)
        k = idx ÷ (nboundary + 1)  # Block index
        i = idx % (nboundary + 1)   # Local row index (0 to nboundary)
        offset_row = k * nboundary
        offset_nnz = k * nnz_b
        # +1 for Julia 1-based indexing
        @inbounds d_rowptr[offset_row + i + 1] = (
            d_rowptrLocal[i + 1] + offset_nnz
        )
    end

    # Each thread handles one colidx entry (offset to avoid overlap)
    idx_col = idx - nb * (nboundary + 1)
    if idx_col >= 0 && idx_col < nb * nnz_b
        k = idx_col ÷ nnz_b
        i = idx_col % nnz_b
        offset_nnz = k * nnz_b
        offset_col = k * numNodes
        @inbounds d_colidx[offset_nnz + i + 1] = (
            d_colidxLocal[i + 1] + offset_col
        )
    end

    return nothing
end

# CUDA kernel to initialize d_utmp from Q1 basis functions
function init_utmp_kernel!(
    d_utmp, d_phiLocal, d_freeToGlobal, nb, nfree, numVectorsToSolve
)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if idx <= nb * nfree
        # Decode: which element and which local DOF
        iel = (idx - 1) ÷ nfree + 1
        i_local = (idx - 1) % nfree + 1

        offset = (iel - 1) * nfree
        i_global = d_freeToGlobal[i_local]

        # Copy the 3 basis function values for this interior node
        for ir = 1:numVectorsToSolve
            d_utmp[offset + i_local, ir] = d_phiLocal[i_global, ir]
        end
    end
    return nothing
end

# CUDA kernel to transfer solution from d_utmp back to d_phi
function transfer_utmp_to_phi_kernel!(
    d_phi, d_utmp, d_freeToGlobal, nb, nfree, numNodes, numVectorsToSolve
)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if idx <= nb * nfree
        # Decode: which block and which local DOF
        iel = (idx - 1) ÷ nfree + 1
        i_local = (idx - 1) % nfree + 1

        # Offset into d_utmp for this block
        utmp_offset = (iel - 1) * nfree

        # Offset into d_phi for this block
        phi_offset = (iel - 1) * numNodes

        # Global node ID within the local element
        i_global = d_freeToGlobal[i_local]

        # Copy the solved values (3 basis functions) from d_utmp to d_phi
        for ir = 1:numVectorsToSolve
            d_phi[phi_offset + i_global, ir] = (
                d_utmp[utmp_offset + i_local, ir]
            )
        end

        # Compute 4th basis function from partition of unity:
        # phi_4 = 1 - phi_1 - phi_2 - phi_3
        d_phi[phi_offset + i_global, 4] = (
            1.0 - d_phi[phi_offset + i_global, 1] -
            d_phi[phi_offset + i_global, 2] -
            d_phi[phi_offset + i_global, 3]
        )
    end
    return nothing
end

# CUDA kernel to reconstruct fine-scale solution from coarse solution
# and basis functions
function reconstruct_fine_solution_kernel!(
    d_fine_x, d_fine_y, d_fine_u, d_coarse_solution,
    d_phi, d_cell_to_node, nx_coarse, ny_coarse, ratio, numNodes
)
    # Global thread index for fine mesh nodes
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    nx_fine = nx_coarse * ratio
    ny_fine = ny_coarse * ratio
    nfine = (nx_fine + 1) * (ny_fine + 1)

    if idx <= nfine
        # Convert linear index to 2D fine mesh coordinates
        idx_0 = idx - 1
        ix_fine = idx_0 % (nx_fine + 1)
        iy_fine = idx_0 ÷ (nx_fine + 1)

        # Store fine mesh coordinates
        @inbounds d_fine_x[idx] = Float64(ix_fine) / nx_fine
        @inbounds d_fine_y[idx] = Float64(iy_fine) / ny_fine

        # Determine which coarse element this fine node belongs to
        ix_coarse = ix_fine ÷ ratio
        iy_coarse = iy_fine ÷ ratio

        # Clamp to coarse mesh bounds (boundary nodes)
        ix_coarse = min(ix_coarse, nx_coarse - 1)
        iy_coarse = min(iy_coarse, ny_coarse - 1)

        # Coarse element index (1-based)
        iel = ix_coarse + iy_coarse * nx_coarse + 1

        # Local fine node coordinates within the coarse element
        ix_local = ix_fine - ix_coarse * ratio
        iy_local = iy_fine - iy_coarse * ratio
        local_node = ix_local + iy_local * (ratio + 1) + 1

        # Get coarse element nodes (Q1: 4 nodes)
        # Offset into d_phi for this element
        offset_phi = (iel - 1) * numNodes

        # Interpolate: u_fine = sum_i phi[local_node, i] * u_coarse[nodes[i]]
        u_val = 0.0
        for i = 1:4
            # Get global coarse node index
            @inbounds node_global = d_cell_to_node[iel, i]
            @inbounds u_coarse_i = d_coarse_solution[node_global]
            @inbounds phi_val = d_phi[offset_phi + local_node, i]
            u_val += phi_val * u_coarse_i
        end

        @inbounds d_fine_u[idx] = u_val
    end

    return nothing
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
            gi = nodes[i]
            fe_val = @inbounds d_fe_coarse[iel, i]
            @inbounds d_rhs[gi] += fe_val
        end

        # Scatter stiffness: mat_values[k] += Ke_coarse[iel, i, j]
        for i = 1:4
            gi = nodes[i]
            for j = 1:4
                gj = nodes[j]

                # Find position in sparse matrix (CSC format)
                # K[gi, gj] is stored at column gj, row gi
                col_start = @inbounds d_col_ptr[gj]
                col_end = @inbounds d_col_ptr[gj + 1] - 1

                # Binary search for row gi in column gj
                k = binary_search_sparse(d_row_idx, gi, col_start, col_end)

                if k != -1
                    # Ke_coarse is stored as nb × 16 (flattened 4×4 mat)
                    ke_idx = (i - 1) * 4 + j
                    ke_val = @inbounds d_Ke_coarse[iel, ke_idx]
                    @inbounds d_mat_values[k] += ke_val
                end
            end
        end
    end

    return nothing
end

# CUDA kernel to compute coarse element RHS and stiffness:
# fe_coarse[iel, ir] = phi[iel, :, ir]' * rhs_fine[iel, :]
# Ke_coarse[iel, ir, jr] = phi[iel, :, ir]' * K_b[iel, :, :] * phi[iel, :, jr]
# Uses shared memory reduction for efficiency and to minimize atomic
# contention. Optimized version caches phi values for boundary nodes
# in shared memory to reduce global memory traffic.
function compute_fe_Ke_coarse_kernel!(
    d_fe_coarse, d_Ke_coarse, d_phi, d_rhs_fine,
    d_valK_b, d_rowptr_b, d_colidx_b, d_boundaryToGlobal,
    nb, numNodes, numVectors, nboundary, nnz_b
)
    # Each block handles one coarse element
    iel = blockIdx().x
    tid = threadIdx().x

    if iel > nb
        return
    end

    # Offsets for this element's data
    offset_nodes = (iel - 1) * numNodes

    # Shared memory layout (all offsets in Float64 elements):
    # - sdata_fe: blockDim().x * numVectors for fe_coarse reduction
    # - sdata_Ke: blockDim().x * 16 for Ke_coarse local accumulation
    # - sdata_phi_boundary: nboundary * numVectors for boundary phi cache
    sdata_fe = @cuDynamicSharedMem(Float64, blockDim().x * numVectors)
    offset_Ke = blockDim().x * numVectors
    sdata_Ke = @cuDynamicSharedMem(
        Float64, blockDim().x * 16,
        offset_Ke * sizeof(Float64)
    )
    offset_phi = offset_Ke + blockDim().x * 16
    sdata_phi_boundary = @cuDynamicSharedMem(
        Float64, nboundary * numVectors,
        offset_phi * sizeof(Float64)
    )

    # Initialize Ke_coarse local accumulator
    for i = 1:16
        @inbounds sdata_Ke[(tid - 1) * 16 + i] = 0.0
    end

    # ====================================================================
    # Load phi values for boundary nodes into shared memory
    # ====================================================================
    # Layout: sdata_phi_boundary[(ib-1)*numVectors + ir] = phi[boundary_node_ib, ir]
    # This avoids repeated global memory reads in the Ke computation
    total_phi_entries = nboundary * numVectors
    for idx = tid:blockDim().x:total_phi_entries
        ib = (idx - 1) ÷ numVectors + 1   # Boundary node index (1 to nboundary)
        ir = (idx - 1) % numVectors + 1   # Basis function index (1 to numVectors)
        @inbounds iGlobal = d_boundaryToGlobal[ib]
        @inbounds sdata_phi_boundary[idx] = d_phi[offset_nodes + iGlobal, ir]
    end

    sync_threads()

    # ====================================================================
    # Part 1: Compute fe_coarse using reduction
    # ====================================================================

    # Each thread computes partial sums for all 4 basis functions
    for ir = 1:numVectors
        local_sum = 0.0
        # Stride loop: each thread processes multiple nodes if
        # numNodes > blockDim().x
        i = tid
        while i <= numNodes
            @inbounds local_sum += (
                d_phi[offset_nodes + i, ir] *
                d_rhs_fine[offset_nodes + i]
            )
            i += blockDim().x
        end
        @inbounds sdata_fe[(tid - 1) * numVectors + ir] = local_sum
    end

    sync_threads()

    # Reduction in shared memory for fe_coarse
    s = blockDim().x ÷ 2
    while s > 0
        if tid <= s
            for ir = 1:numVectors
                idx = (tid - 1) * numVectors + ir
                idx_s = (tid + s - 1) * numVectors + ir
                @inbounds sdata_fe[idx] += sdata_fe[idx_s]
            end
        end
        sync_threads()
        s = s ÷ 2
    end

    # Thread 1 writes fe_coarse result
    if tid == 1
        for ir = 1:numVectors
            @inbounds d_fe_coarse[iel, ir] = sdata_fe[ir]
        end
    end

    sync_threads()

    # ====================================================================
    # Part 2: Compute Ke_coarse = phi^T * K_b * phi (local accumulation)
    # ====================================================================

    # Parallelize over boundary rows: each thread handles some rows
    # Accumulate locally in shared memory to avoid excessive atomics

    # Block offset for global boundary indexing
    block_offset_boundary = (iel - 1) * nboundary

    for iBoundary_local = tid:blockDim().x:nboundary
        if iBoundary_local <= nboundary
            # Global boundary index for rowptr lookup
            iBoundary_global = block_offset_boundary + iBoundary_local

            # Row range in K_b for this boundary node (using global index)
            @inbounds row_start = d_rowptr_b[iBoundary_global]
            @inbounds row_end = d_rowptr_b[iBoundary_global + 1] - 1

            # Iterate over non-zeros in this row
            for k_pos = row_start:row_end
                # k_pos is already global, no offset needed
                @inbounds jGlobal_offset = d_colidx_b[k_pos]
                @inbounds k_val = d_valK_b[k_pos]

                # Convert global offset to local node index
                jGlobal = jGlobal_offset - offset_nodes

                # Accumulate contribution to Ke_coarse for all basis
                # function pairs. Store in thread-local shared memory
                # (NO atomics here!)
                for ir = 1:numVectors
                    # Read phi_i from shared memory (cached boundary values)
                    phi_idx_i = (iBoundary_local - 1) * numVectors + ir
                    @inbounds phi_i = sdata_phi_boundary[phi_idx_i]
                    phi_i_k = phi_i * k_val

                    for jr = 1:numVectors
                        # Read phi_j from global memory (jGlobal can be any node)
                        @inbounds phi_j = d_phi[offset_nodes + jGlobal, jr]
                        contribution = phi_i_k * phi_j

                        # Linear index into local Ke_coarse (4x4 = 16 elems)
                        ke_idx = (ir - 1) * numVectors + jr
                        # Accumulate in shared memory (thread-private)
                        @inbounds sdata_Ke[(tid - 1) * 16 + ke_idx] += (
                            contribution
                        )
                    end
                end
            end
        end
    end

    sync_threads()

    # Final reduction: sum all thread-local Ke_coarse contributions
    # Only use atomics at the very end (16 atomics per thread instead
    # of thousands!)
    for i = 1:16
        local_val = sdata_Ke[(tid - 1) * 16 + i]
        if local_val != 0.0
            CUDA.@atomic d_Ke_coarse[iel, i] += local_val
        end
    end

    return nothing
end

# ========================================================================
# Main Program
# ========================================================================

function main()
    println("="^70)
    println("Julia 2D MFEM Assembly (Cuda with Coloring)")
    println("="^70)
    println()

    # ====================================================================
    # Mesh generation
    # ====================================================================

    nx, ny = 32, 32  # Coarse mesh (MFEM)
    ratio = 32        # Fine grid refinement ratio per coarse element

    # Validate ratio is a power of 2 (required for exact floating-point arithmetic)
    if !ispow2(ratio)
        error("ratio must be a power of 2 (got $ratio). Use 16, 32, 64, 128, etc.")
    end
    log2_ratio = trailing_zeros(ratio)  # Fast log2 for power of 2

    println("Generating coarse mesh: $nx x $ny MFEM elements")
    println(
        "MFEM ratio: $ratio (each coarse element has " *
        "$(ratio)x$(ratio) fine elements)"
    )
    println("Effective fine resolution: $(nx*ratio) x $(ny*ratio)")
    println()

    t0 = time()
    mesh = generate_mesh(nx, ny)
    mesh_time = time() - t0

    println("  Number of coarse elements: ", size(mesh.cell_to_node, 1))
    println("  Number of coarse nodes:    ", length(mesh.vertex_x))
    println(
        "  Mesh generation:           ",
        @sprintf("%.2f ms", mesh_time * 1000)
    )
    println()

    # ====================================================================
    # Build connectivity and coloring
    # ====================================================================

    # Warm-up run to eliminate JIT compilation overhead
    println("Running connectivity JIT warm-up...")
    t_warmup = time()
    _, _, _ = build_mesh_connectivity(mesh)
    warmup_time = time() - t_warmup
    println("  Warm-up time: ", @sprintf("%.2f ms", warmup_time * 1000))
    println()

    println("Building mesh connectivity...")
    t0 = time()
    n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
    conn_time = time() - t0

    nnz = length(n2n_col_idx)
    println(
        "  Matrix size: ", length(mesh.vertex_x), " x ",
        length(mesh.vertex_x)
    )
    println("  Non-zeros:   ", nnz)
    println("  Number of colors: ", length(e2e_colors))
    for ic = 1:length(e2e_colors)
        println("    Color $ic: ", length(e2e_colors[ic]), " elements")
    end
    println(
        "  Coarse Connectivity time: ",
        @sprintf("%.2f ms", conn_time * 1000)
    )
    println()

    # Build the connectivity for one coarse element

    println("Building local mesh connectivity...")
    t0 = time()
    mesh_local = generate_mesh(ratio, ratio)
    n2n_row_ptr_local, n2n_col_idx_local, _ = build_mesh_connectivity(
        mesh_local
    )
    local_conn_time = time() - t0

    nnz = length(n2n_col_idx_local)
    println(
        "  Local Matrix size: ", length(mesh_local.vertex_x), " x ",
        length(mesh_local.vertex_x)
    )
    println("  Non-zeros:   ", nnz)
    println(
        "  Local Connectivity time: ",
        @sprintf("%.2f ms", local_conn_time * 1000)
    )
    println()

    # ====================================================================
    # JIT Warm-up for GPU Kernels
    # ====================================================================

    println("Running GPU kernel JIT warm-up...")
    t_warmup_gpu = time()

    # Create small dummy problem for warm-up (1 element, ratio=4)
    nb_warmup = 1
    ratio_warmup = 4
    numNodes_warmup = (ratio_warmup + 1)^2
    workspace_warmup = MFEWorkspace{Float64, 2, 4}(ratio_warmup)
    nfree_w, nboundary_w, nnz_ii_w, nnz_b_w = build_mfem_local_sparsity!(
        workspace_warmup, ratio_warmup, ratio_warmup
    )

    # Warm-up arrays (tiny sizes to minimize overhead)
    d_phi_w = CUDA.zeros(Float64, nb_warmup * numNodes_warmup, 4)
    d_valK_ii_w = CUDA.zeros(Float64, nb_warmup * nnz_ii_w)
    d_colptr_ii_w = CUDA.zeros(Int32, nb_warmup * nfree_w + 1)
    d_rowidx_ii_w = CUDA.zeros(Int32, nb_warmup * nnz_ii_w)
    d_btmp_w = CUDA.zeros(Float64, nb_warmup * nfree_w, 3)
    d_rhs_fine_w = CUDA.zeros(Float64, nb_warmup * numNodes_warmup)
    d_valK_b_w = CUDA.zeros(Float64, nb_warmup * nnz_b_w)
    d_rowptr_b_w = CUDA.zeros(Int32, nb_warmup * nboundary_w + 1)
    d_colidx_b_w = CUDA.zeros(Int32, nb_warmup * nnz_b_w)
    d_globalToFree_w = cu(workspace_warmup.globalToFree)
    d_globalToBoundary_w = cu(workspace_warmup.globalToBoundary)
    d_freeToGlobal_w = cu(workspace_warmup.freeToGlobal)
    d_boundaryToGlobal_w = cu(workspace_warmup.boundaryToGlobal)

    # Initialize d_colptr and d_rowidx for warm-up (valid sparsity)
    d_colptrLocal_w = cu(workspace_warmup.colptr_ii)
    d_rowidxLocal_w = cu(workspace_warmup.rowidx_ii)
    d_rowptrLocal_b_w = cu(workspace_warmup.matRowPtr_b)
    d_colidxLocal_b_w = cu(workspace_warmup.matColIdx_b)

    # Warm-up kernel 1: build_block_diagonal_sparsity
    total_work_w = nb_warmup * (nfree_w + 1) + nb_warmup * nnz_ii_w
    @cuda threads=64 blocks=cld(total_work_w, 64) (
        build_block_diagonal_sparsity_kernel!(
            d_colptr_ii_w, d_rowidx_ii_w, d_colptrLocal_w, d_rowidxLocal_w,
            nb_warmup, nfree_w, nnz_ii_w
        )
    )
    CUDA.synchronize()

    # Warm-up kernel 2: build_block_diagonal_csr_sparsity
    total_work_b_w = nb_warmup * (nboundary_w + 1) + nb_warmup * nnz_b_w
    @cuda threads=64 blocks=cld(total_work_b_w, 64) (
        build_block_diagonal_csr_sparsity_kernel!(
            d_rowptr_b_w, d_colidx_b_w, d_rowptrLocal_b_w,
            d_colidxLocal_b_w, nb_warmup, nboundary_w, nnz_b_w,
            numNodes_warmup
        )
    )
    CUDA.synchronize()

    # Warm-up kernel 3: init_utmp
    d_utmp_w = CUDA.zeros(Float64, nb_warmup * nfree_w, 3)
    d_phiLocal_w = cu(workspace_warmup.phi)
    @cuda threads=64 blocks=cld(nb_warmup * nfree_w, 64) (
        init_utmp_kernel!(
            d_utmp_w, d_phiLocal_w, d_freeToGlobal_w,
            nb_warmup, nfree_w, 3
        )
    )
    CUDA.synchronize()

    # Warm-up kernel 4: assemble_Kii (main assembly kernel)
    log2_ratio_warmup = trailing_zeros(ratio_warmup)
    @cuda threads=64 blocks=1 assemble_Kii_kernel!(
        d_valK_ii_w, d_btmp_w, d_rhs_fine_w, d_phi_w, d_colptr_ii_w,
        d_rowidx_ii_w, d_globalToFree_w, d_valK_b_w, d_rowptr_b_w,
        d_colidx_b_w, d_globalToBoundary_w, 1, 1, ratio_warmup, log2_ratio_warmup, nfree_w,
        nnz_ii_w, nboundary_w, nnz_b_w, 3, numNodes_warmup
    )
    CUDA.synchronize()

    # Warm-up kernel 5: block_cg solver (critical!)
    # Need a valid sparse matrix structure - build on CPU then upload
    valK_ii_w_cpu = ones(Float64, nfree_w > 0 ? length(workspace_warmup.nzval_ii) : 1)
    for i = 1:nfree_w
        col_start = workspace_warmup.colptr_ii[i]
        col_end = workspace_warmup.colptr_ii[i+1] - 1
        for k = col_start:col_end
            if workspace_warmup.rowidx_ii[k] == i
                valK_ii_w_cpu[k] = 10.0  # Diagonal - make dominant
            end
        end
    end
    copyto!(d_valK_ii_w, valK_ii_w_cpu[1:min(length(valK_ii_w_cpu), nb_warmup * nnz_ii_w)])
    K_warmup = CuSparseMatrixCSC(
        d_colptr_ii_w, d_rowidx_ii_w, d_valK_ii_w, (nfree_w, nfree_w)
    )
    d_btmp_w .= 1.0
    # Warm-up extract_diagonal_kernel
    d_diag_w = CUDA.zeros(Float64, nfree_w)
    @cuda threads=64 blocks=cld(nfree_w, 64) extract_diagonal_kernel!(
        d_diag_w, d_colptr_ii_w, d_rowidx_ii_w, d_valK_ii_w, nfree_w
    )
    CUDA.synchronize()
    # Build Jacobi preconditioner for warm-up
    d_precond_w = map(d -> d ≠ 0 ? 1 / abs(d) : 1.0, d_diag_w)
    precond_w = Diagonal(d_precond_w)
    # Use regular CG for each RHS column (warm-up)
    for col = 1:size(d_btmp_w, 2)
        x_col = view(d_utmp_w, :, col)
        b_col = view(d_btmp_w, :, col)
        cg(K_warmup, b_col; M=precond_w, atol=1e-6, rtol=1e-6, itmax=5, verbose=0)
    end
    CUDA.synchronize()

    # Warm-up kernel 6: transfer_utmp_to_phi
    @cuda threads=64 blocks=cld(nb_warmup * nfree_w, 64) (
        transfer_utmp_to_phi_kernel!(
            d_phi_w, d_utmp_w, d_freeToGlobal_w,
            nb_warmup, nfree_w, numNodes_warmup, 3
        )
    )
    CUDA.synchronize()

    # Warm-up kernel 7: compute_fe_Ke_coarse
    d_fe_coarse_w = CUDA.zeros(Float64, nb_warmup, 4)
    d_Ke_coarse_w = CUDA.zeros(Float64, nb_warmup, 16)
    # Shared memory: fe (32*4) + Ke (32*16) + phi_boundary (nboundary_w*4)
    shared_mem_w = (32 * (4 + 16) + nboundary_w * 4) * sizeof(Float64)
    @cuda threads=32 blocks=1 shmem=shared_mem_w (
        compute_fe_Ke_coarse_kernel!(
            d_fe_coarse_w, d_Ke_coarse_w, d_phi_w, d_rhs_fine_w,
            d_valK_b_w, d_rowptr_b_w, d_colidx_b_w,
            d_boundaryToGlobal_w, nb_warmup, numNodes_warmup, 4,
            nboundary_w, nnz_b_w
        )
    )
    CUDA.synchronize()

    # Warm-up kernel 8: reconstruct_fine_solution
    d_coarse_sol_w = CUDA.zeros(Float64, 4)
    d_cell_to_node_w = cu(reshape([1, 2, 4, 3], 1, 4))
    d_fine_x_w = CUDA.zeros(Float64, 25)
    d_fine_y_w = CUDA.zeros(Float64, 25)
    d_fine_u_w = CUDA.zeros(Float64, 25)
    @cuda threads=64 blocks=1 reconstruct_fine_solution_kernel!(
        d_fine_x_w, d_fine_y_w, d_fine_u_w, d_coarse_sol_w,
        d_phi_w, d_cell_to_node_w, 1, 1, ratio_warmup, numNodes_warmup
    )
    CUDA.synchronize()

    # Warm-up kernel 9: assemble_coarse_color (new!)
    d_mat_values_w = CUDA.zeros(Float64, 16)
    d_rhs_w = CUDA.zeros(Float64, 4)
    d_col_ptr_w = cu(Int32[1, 5, 9, 13, 17])  # 4 nodes, 4 nnz per node
    d_row_idx_w = cu(Int32[
        1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4
    ])
    d_elements_w = cu(Int32[1])  # Single element
    @cuda threads=64 blocks=1 assemble_coarse_color_kernel!(
        d_mat_values_w, d_rhs_w, d_fe_coarse_w, d_Ke_coarse_w,
        d_cell_to_node_w, d_col_ptr_w, d_row_idx_w,
        d_elements_w, 1
    )
    CUDA.synchronize()

    warmup_time_gpu = time() - t_warmup_gpu
    println(
        "  GPU warm-up time: ",
        @sprintf("%.2f ms", warmup_time_gpu * 1000)
    )
    println("  (All 9 CUDA kernels and block CG solver now compiled)")
    println()

    workspace = MFEWorkspace{Float64, 2, 4}(ratio)

    # Build DOF mappings and sparsity patterns using reusable kernel
    # from FEMBase
    nfree, nboundary, nnz_ii, nnz_b = build_mfem_local_sparsity!(
        workspace, ratio, ratio
    )

    # Extract arrays from workspace for later use
    globalToFree = workspace.globalToFree
    freeToGlobal = workspace.freeToGlobal

    # Fill phi with Q1 basis function values on the fine grid
    # phi[nodeID, basis] = value of basis function 'basis' at node 'nodeID'
    phi = workspace.phi
    fill!(phi, 0.0)

    for iy = 0:ratio
        for ix = 0:ratio
            # Node ID in the fine grid (1-based indexing)
            nodeID = ix + iy * (ratio + 1) + 1

            # Normalized coordinates within the coarse element [0,1] × [0,1]
            # Use ldexp for exact division by power of 2
            xi = ldexp(Float64(ix), -log2_ratio)
            eta = ldexp(Float64(iy), -log2_ratio)

            # Evaluate Q1 basis functions at (xi, eta)
            # Q1 nodes: [1] bottom-left, [2] bottom-right,
            #           [3] top-right, [4] top-left
            @inbounds phi[nodeID, 1] = (1.0 - xi) * (1.0 - eta)  # Bottom-left
            @inbounds phi[nodeID, 2] = xi * (1.0 - eta)          # Bottom-right
            @inbounds phi[nodeID, 3] = xi * eta                  # Top-right
            @inbounds phi[nodeID, 4] = (1.0 - xi) * eta          # Top-left
        end
    end

    nb = nx * ny  # Number of blocks
    numNodes = (ratio + 1) * (ratio + 1)  # Fine nodes per coarse element

    d_phiLocal = cu(phi)
    d_phi = repeat(d_phiLocal, nb, 1)

    # Build column pointer for block diagonal sparse matrix
    # Each of the nx*ny blocks has the same sparsity pattern as K_ii

    d_colptr_ii = CUDA.zeros(Int32, nb * nfree + 1)
    d_rowidx_ii = CUDA.zeros(Int32, nb * nnz_ii)

    # Build block diagonal sparsity pattern - on device
    d_colptrLocal = cu(workspace.colptr_ii)
    d_rowidxLocal = cu(workspace.rowidx_ii)

    # Calculate total work: colptr entries + rowidx entries
    threads_per_block = 256
    total_work_ii = nb * (nfree + 1) + nb * nnz_ii
    num_blocks = cld(total_work_ii, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks (
        build_block_diagonal_sparsity_kernel!(
            d_colptr_ii, d_rowidx_ii, d_colptrLocal, d_rowidxLocal,
            nb, nfree, nnz_ii
        )
    )
    CUDA.synchronize()

    d_valK_ii = CUDA.zeros(Float64, nb * nnz_ii)

    # Build block diagonal CSR sparsity pattern for K_b
    # (boundary-all coupling)
    # K_b has nboundary rows and numNodes columns per element
    d_rowptr_b = CUDA.zeros(Int32, nb * nboundary + 1)
    d_colidx_b = CUDA.zeros(Int32, nb * nnz_b)

    # Upload local K_b sparsity pattern from workspace
    d_rowptrLocal_b = cu(workspace.matRowPtr_b)
    d_colidxLocal_b = cu(workspace.matColIdx_b)

    # Calculate total work: rowptr entries + colidx entries
    threads_per_block = 256
    total_work_b = nb * (nboundary + 1) + nb * nnz_b
    num_blocks = cld(total_work_b, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks (
        build_block_diagonal_csr_sparsity_kernel!(
            d_rowptr_b, d_colidx_b, d_rowptrLocal_b, d_colidxLocal_b,
            nb, nboundary, nnz_b, numNodes
        )
    )
    CUDA.synchronize()

    d_valK_b = CUDA.zeros(Float64, nb * nnz_b)

    # Allocate btmp (RHS for static condensation):
    # (nb*nfree) × numVectorsToSolve
    numVectorsToSolve = 3
    d_btmp = CUDA.zeros(Float64, nb * nfree, numVectorsToSolve)

    # Upload globalToFree and globalToBoundary mappings to device
    d_globalToFree = cu(workspace.globalToFree)
    d_globalToBoundary = cu(workspace.globalToBoundary)

    # Allocate fine-grid RHS for all nodes (interior + boundary)
    # in all coarse elements
    d_rhs_fine = CUDA.zeros(Float64, nb * numNodes)

    # Launch kernel to assemble K_ii, K_b, btmp, and rhs_fine
    println("GPU fine-grid assembly...")
    t0 = time()
    threads_per_block = 256
    num_blocks = cld(nb, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks assemble_Kii_kernel!(
        d_valK_ii, d_btmp, d_rhs_fine, d_phi, d_colptr_ii, d_rowidx_ii,
        d_globalToFree, d_valK_b, d_rowptr_b, d_colidx_b,
        d_globalToBoundary, nx, ny, ratio, log2_ratio, nfree, nnz_ii, nboundary,
        nnz_b, numVectorsToSolve, numNodes
    )
    CUDA.synchronize()
    gpu_fine_assembly_time = time() - t0
    println(
        "  GPU fine assembly time: ",
        @sprintf("%.2f ms", gpu_fine_assembly_time * 1000)
    )
    println()

    # Create the solution vector and initial guess
    d_utmp = CUDA.zeros(Float64, nb * nfree, numVectorsToSolve)
    d_freeToGlobal = cu(workspace.freeToGlobal)

    # Build initial guess from Q1 basis functions (partition of unity)
    # - on device
    threads_per_block = 256
    num_blocks = cld(nb * nfree, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks init_utmp_kernel!(
        d_utmp, d_phiLocal, d_freeToGlobal, nb, nfree, numVectorsToSolve
    )
    CUDA.synchronize()

    # Build the block diagonal sparse matrix K_ii on device
    n_total = nb * nfree
    K_ii_gpu = CuSparseMatrixCSC(
        d_colptr_ii, d_rowidx_ii, d_valK_ii, (n_total, n_total)
    )

    # Build diagonal Jacobi preconditioner for K_ii
    # P⁻¹ = diag(1/|K[i,i]|)
    println("Building diagonal Jacobi preconditioner...")
    t0_prec = time()
    # Extract diagonal from sparse matrix on GPU
    d_diag = CUDA.zeros(Float64, n_total)
    @cuda threads=256 blocks=cld(n_total, 256) extract_diagonal_kernel!(
        d_diag, d_colptr_ii, d_rowidx_ii, d_valK_ii, n_total
    )
    CUDA.synchronize()
    # Build preconditioner: P⁻¹[i] = 1/|diag[i]| (or 1 if zero)
    d_precond = map(d -> d ≠ 0 ? 1 / abs(d) : 1.0, d_diag)
    precond = Diagonal(d_precond)
    prec_time = time() - t0_prec
    println("  Preconditioner build time: ", @sprintf("%.2f ms", prec_time * 1000))

    # Solve K_ii * d_utmp = d_btmp using CG for each RHS column
    # Use CgSolver workspace for memory reuse and proper warm-start
    println("Solving CG on GPU (3 RHS columns)...")
    t0 = time()
    total_iters = 0
    all_converged = true

    # Create CG solver workspace once (reused for all RHS columns)
    b_col = view(d_btmp, :, 1)
    solver = CgWorkspace(K_ii_gpu, b_col)

    for col = 1:size(d_btmp, 2)
        b_col = view(d_btmp, :, col)
        x0_col = view(d_utmp, :, col)

        # Check if initial guess is close to solution
        r0 = b_col - K_ii_gpu * x0_col
        r0_norm = norm(r0)
        b_norm = norm(b_col)
        rel_residual = r0_norm / max(b_norm, 1e-16)
        println("  RHS $col: ||r0|| = ", @sprintf("%.2e", r0_norm),
                ", ||b|| = ", @sprintf("%.2e", b_norm),
                ", relative = ", @sprintf("%.2e", rel_residual))

        # Skip CG if initial guess is already converged
        if rel_residual < 1e-12
            # Already converged - no iterations needed
            # Solution is already in x0_col (d_utmp), no copy needed
            total_iters += 0
        else
            # Need to solve - use cg! with warm-start and block Jacobi preconditioner
            solver = cg!(solver, K_ii_gpu, b_col, x0_col;
                M=precond, atol=1e-24, rtol=1e-12, itmax=1000, verbose=0)
            copyto!(view(d_utmp, :, col), solver.x)
            total_iters += solver.stats.niter
            all_converged = all_converged && solver.stats.solved
        end
    end
    gpu_cg_time = time() - t0

    println("  CG converged: ", all_converged)
    println("  Total iterations: ", total_iters)
    println("  GPU CG time: ", @sprintf("%.2f ms", gpu_cg_time * 1000))
    println()

    # Transfer solution from d_utmp back to d_phi
    threads_per_block = 256
    num_blocks = cld(nb * nfree, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks (
        transfer_utmp_to_phi_kernel!(
            d_phi, d_utmp, d_freeToGlobal, nb, nfree, numNodes,
            numVectorsToSolve
        )
    )
    CUDA.synchronize()

    # Compute coarse element RHS and stiffness matrices
    # fe_coarse[iel, ir] = phi[iel, :, ir]' * rhs_fine[iel, :]
    # Ke_coarse[iel, ir*4+jr] = phi[iel, :, ir]' * K_b * phi[iel, :, jr]
    println("GPU coarse element assembly...")
    t0 = time()
    numVectors = 4
    d_fe_coarse = CUDA.zeros(Float64, nb, numVectors)
    d_Ke_coarse = CUDA.zeros(Float64, nb, numVectors * numVectors)

    # Upload boundaryToGlobal mapping to device
    d_boundaryToGlobal = cu(workspace.boundaryToGlobal)

    threads_per_block = 128  # Use 128 threads for better occupancy
    num_blocks = nb  # One block per coarse element
    # Shared memory layout:
    # - fe_coarse reduction: threads * numVectors
    # - Ke_coarse local accumulation: threads * 16
    # - phi boundary cache: nboundary * numVectors (NEW)
    shared_mem_size = (
        threads_per_block * (numVectors + 16) + nboundary * numVectors
    ) * sizeof(Float64)

    @cuda threads=threads_per_block blocks=num_blocks shmem=shared_mem_size (
        compute_fe_Ke_coarse_kernel!(
            d_fe_coarse, d_Ke_coarse, d_phi, d_rhs_fine, d_valK_b,
            d_rowptr_b, d_colidx_b, d_boundaryToGlobal, nb, numNodes,
            numVectors, nboundary, nnz_b
        )
    )
    CUDA.synchronize()
    gpu_coarse_assembly_time = time() - t0
    println(
        "  GPU coarse assembly time: ",
        @sprintf("%.2f ms", gpu_coarse_assembly_time * 1000)
    )
    println()

    # ====================================================================
    # GPU Assembly for the coarse problem
    # ====================================================================

    println("="^70)
    println("Starting MFEM Coarse Assembly on GPU")
    println("="^70)

    println("Assembling coarse system on GPU using coloring...")
    t0 = time()

    # Allocate global matrix and RHS on GPU
    d_mat_values = CUDA.zeros(Float64, length(n2n_col_idx))
    d_rhs = CUDA.zeros(Float64, length(mesh.vertex_x))

    # Upload sparse matrix structure and connectivity to GPU
    # Note: despite name, this is column pointer (CSC format)
    d_col_ptr = cu(n2n_row_ptr)
    # Note: despite name, this is row index (CSC format)
    d_row_idx = cu(n2n_col_idx)
    d_cell_to_node = cu(mesh.cell_to_node)

    # Assemble by color to avoid race conditions
    for (ic, elements) in enumerate(e2e_colors)
        # Upload elements for this color
        d_elements = cu(Int32.(elements))
        num_elements = length(elements)

        # Launch kernel: one thread per element in this color
        threads_per_block = 256
        num_blocks = cld(num_elements, threads_per_block)

        @cuda threads=threads_per_block blocks=num_blocks (
            assemble_coarse_color_kernel!(
                d_mat_values, d_rhs, d_fe_coarse, d_Ke_coarse,
                d_cell_to_node, d_col_ptr, d_row_idx,
                d_elements, num_elements
            )
        )
    end
    CUDA.synchronize()

    assembly_time = time() - t0
    println("  Number of colors: ", length(e2e_colors))
    println(
        "  GPU assembly time: ",
        @sprintf("%.2f ms", assembly_time * 1000)
    )

    # Transfer assembled matrix and RHS to CPU
    mat_values = Array(d_mat_values)
    rhs = Array(d_rhs)
    println("  Transfer time included in assembly time")
    println()

    # ====================================================================
    # Apply boundary conditions
    # ====================================================================

    println("Applying boundary conditions...")
    t0 = time()
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs,
        free_to_global = apply_boundary_conditions(
            n2n_row_ptr, n2n_col_idx, mat_values, rhs, mesh.boundary_nodes
        )
    bc_time = time() - t0

    nfree = length(free_to_global)
    println("  Total DOFs:    ", length(mesh.vertex_x))
    println("  Boundary DOFs: ", length(mesh.boundary_nodes))
    println("  Free DOFs:     ", nfree)
    println("  Reduced nnz:   ", length(reduced_col_idx))
    println("  BC time:       ", @sprintf("%.2f ms", bc_time * 1000))
    println()

    # ====================================================================
    # Solve
    # ====================================================================

    println("Solving linear system...")

    # Build sparse matrix
    K = SparseMatrixCSC(
        nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values
    )

    t0 = time()
    sol_free = K \ reduced_rhs
    solve_time = time() - t0

    println("  Solve time: ", @sprintf("%.2f ms", solve_time * 1000))
    println()

    # ====================================================================
    # Reconstruct Fine-Scale Solution (on GPU)
    # ====================================================================

    println("Reconstructing fine-scale solution on GPU...")

    # Expand coarse solution to full DOFs (BCs = 0) and upload to GPU
    coarse_solution = zeros(length(mesh.vertex_x))
    for i = 1:nfree
        @inbounds coarse_solution[free_to_global[i]] = sol_free[i]
    end
    d_coarse_solution = cu(coarse_solution)

    # Upload cell_to_node connectivity to GPU
    d_cell_to_node = cu(mesh.cell_to_node)

    # Allocate fine mesh arrays on GPU
    nx_fine = nx * ratio
    ny_fine = ny * ratio
    nfine = (nx_fine + 1) * (ny_fine + 1)
    d_fine_x = CUDA.zeros(Float64, nfine)
    d_fine_y = CUDA.zeros(Float64, nfine)
    d_fine_u = CUDA.zeros(Float64, nfine)

    # Launch reconstruction kernel
    t0 = time()
    threads_per_block = 256
    num_blocks = cld(nfine, threads_per_block)
    @cuda threads=threads_per_block blocks=num_blocks (
        reconstruct_fine_solution_kernel!(
            d_fine_x, d_fine_y, d_fine_u, d_coarse_solution,
            d_phi, d_cell_to_node, nx, ny, ratio, numNodes
        )
    )
    CUDA.synchronize()
    reconstruct_time = time() - t0

    # Transfer fine solution from GPU to CPU (only final result)
    fine_x = Array(d_fine_x)
    fine_y = Array(d_fine_y)
    fine_u = Array(d_fine_u)

    println("  Fine mesh size: ", length(fine_u), " nodes")
    println(
        "  GPU reconstruction time: ",
        @sprintf("%.2f ms", reconstruct_time * 1000)
    )
    println()

    # ====================================================================
    # Output
    # ====================================================================

    println("Writing solutions to files...")

    # Write coarse solution
    write_solution("mfem_solution_coarse.txt", mesh, coarse_solution)
    println("  Coarse solution written to: mfem_solution_coarse.txt")

    # Write fine solution
    open("mfem_solution_fine.txt", "w") do io
        println(io, "# x y u")
        for i = 1:length(fine_u)
            @printf(io, "%.6f %.6f %.6e\n", fine_x[i], fine_y[i], fine_u[i])
        end
    end
    println("  Fine solution written to: mfem_solution_fine.txt")
    println()

    # ====================================================================
    # Performance summary
    # ====================================================================

    println("="^70)
    println("Performance Summary (in execution order)")
    println("="^70)
    println(
        "  Mesh generation:       ",
        @sprintf("%8.2f ms", mesh_time * 1000), "  (CPU)"
    )
    println(
        "  Coarse connectivity:   ",
        @sprintf("%8.2f ms", conn_time * 1000), "  (CPU)"
    )
    println(
        "  Local connectivity:    ",
        @sprintf("%8.2f ms", local_conn_time * 1000), "  (CPU)"
    )
    println(
        "  Fine assembly:         ",
        @sprintf("%8.2f ms", gpu_fine_assembly_time * 1000), "  (GPU)"
    )
    println(
        "  CG solve (3 RHS):      ",
        @sprintf("%8.2f ms", gpu_cg_time * 1000), "  (GPU)"
    )
    println(
        "  Coarse element Ke/fe:  ",
        @sprintf("%8.2f ms", gpu_coarse_assembly_time * 1000), "  (GPU)"
    )
    println(
        "  Coarse scatter:        ",
        @sprintf("%8.2f ms", assembly_time * 1000), "  (GPU)"
    )
    println(
        "  Boundary conditions:   ",
        @sprintf("%8.2f ms", bc_time * 1000), "  (CPU)"
    )
    println(
        "  Coarse solve:          ",
        @sprintf("%8.2f ms", solve_time * 1000), "  (CPU)"
    )
    println(
        "  Reconstruction:        ",
        @sprintf("%8.2f ms", reconstruct_time * 1000), "  (GPU)"
    )
    println()
    println("  " * "-"^68)
    cpu_time = (
        mesh_time + conn_time + local_conn_time + bc_time + solve_time
    )
    gpu_time = (
        gpu_fine_assembly_time + gpu_cg_time +
        gpu_coarse_assembly_time + assembly_time + reconstruct_time
    )
    total_time = cpu_time + gpu_time
    println(
        "  CPU Total:             ",
        @sprintf("%8.2f ms", cpu_time * 1000)
    )
    println(
        "  GPU Total:             ",
        @sprintf("%8.2f ms", gpu_time * 1000)
    )
    println(
        "  Total:                 ",
        @sprintf("%8.2f ms", total_time * 1000)
    )
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
