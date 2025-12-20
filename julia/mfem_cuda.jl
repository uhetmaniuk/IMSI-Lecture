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
using Krylov
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


# ============================================================================
# CUDA Kernels
# ============================================================================

# Device functions for material coefficients and forcing term
@inline function ax_func(x, y)
    return 1.0  # Diffusion coefficient in x-direction
end

@inline function ay_func(x, y)
    return 1.0  # Diffusion coefficient in y-direction
end

@inline function f_func(x, y)
    # Manufactured solution: u = sin(π*x) * sin(π*y)
    return 2.0 * π^2 * sin(π * x) * sin(π * y)
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

# CUDA kernel to assemble K_ii, btmp (RHS), and rhs_fine for all coarse elements
function assemble_Kii_kernel!(d_valK_ii, d_btmp, d_rhs_fine, d_phi, d_colptr, d_rowidx, d_globalToFree,
                               d_valK_b, d_rowptr_b, d_colidx_b, d_globalToBoundary,
                               nx, ny, ratio, nfree, nnz_ii, nboundary, nnz_b, numVectorsToSolve, numNodes)
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
        hx = hx_coarse / ratio
        hy = hy_coarse / ratio

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

                # Assemble Q1 element stiffness matrix (2x2 Gauss quadrature)
                # For axis-aligned rectangular elements, Jacobian is constant
                J11 = 0.5 * hx
                J22 = 0.5 * hy
                detJ = J11 * J22
                invJ11 = 1.0 / J11
                invJ22 = 1.0 / J22

                # Gauss quadrature: 2x2 rule
                gp = 1.0 / sqrt(3.0)

                # Initialize element stiffness and RHS
                Ke_fine = ntuple(i -> 0.0, 16)  # 4x4 matrix as tuple
                fe_fine = ntuple(i -> 0.0, 4)   # 4-vector as tuple

                # Loop over 4 Gauss points
                for (qp, (xi, eta)) in enumerate(((-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)))
                    # Q1 shape functions at Gauss point
                    N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
                    N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
                    N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
                    N4 = 0.25 * (1.0 - xi) * (1.0 + eta)

                    # Physical coordinates at Gauss point
                    # Fine element corners: (x1,y1), (x1+hx,y1), (x1+hx,y1+hy), (x1,y1+hy)
                    xq = x1 * N1 + (x1 + hx) * N2 + (x1 + hx) * N3 + x1 * N4
                    yq = y1 * N1 + y1 * N2 + (y1 + hy) * N3 + (y1 + hy) * N4

                    # Evaluate coefficients at Gauss point
                    ax_val = ax_func(xq, yq)
                    ay_val = ay_func(xq, yq)
                    f_val = f_func(xq, yq)

                    # Q1 shape function gradients in reference coords
                    dN_dxi = (-0.25 * (1.0 - eta), 0.25 * (1.0 - eta),
                              0.25 * (1.0 + eta), -0.25 * (1.0 + eta))
                    dN_deta = (-0.25 * (1.0 - xi), -0.25 * (1.0 + xi),
                               0.25 * (1.0 + xi), 0.25 * (1.0 - xi))

                    # Physical gradients
                    dN_dx = ntuple(i -> dN_dxi[i] * invJ11, 4)
                    dN_dy = ntuple(i -> dN_deta[i] * invJ22, 4)

                    # Quadrature weight * detJ (weight = 1 for 2x2 Gauss)
                    w_detJ = detJ

                    # Shape function values at Gauss point
                    N = (N1, N2, N3, N4)

                    # Accumulate stiffness
                    Ke_fine = ntuple(16) do idx
                        i = (idx - 1) ÷ 4 + 1
                        j = (idx - 1) % 4 + 1
                        Ke_fine[idx] + (ax_val * dN_dx[i] * dN_dx[j] +
                                       ay_val * dN_dy[i] * dN_dy[j]) * w_detJ
                    end

                    # Accumulate RHS: fe[i] = ∫ N_i * f dΩ
                    fe_fine = ntuple(4) do i
                        fe_fine[i] + N[i] * f_val * w_detJ
                    end
                end

                # Fine element node list (local numbering within coarse element)
                nodeList = (ix_fine + iy_fine * (ratio + 1) + 1,
                           ix_fine + 1 + iy_fine * (ratio + 1) + 1,
                           ix_fine + 1 + (iy_fine + 1) * (ratio + 1) + 1,
                           ix_fine + (iy_fine + 1) * (ratio + 1) + 1)

                # Scatter RHS to global d_rhs_fine
                offset_rhs = (iel - 1) * numNodes
                for i = 1:4
                    iGlobal = nodeList[i]
                    rhs_idx = offset_rhs + iGlobal
                    CUDA.@atomic d_rhs_fine[rhs_idx] += fe_fine[i]
                end

                # Scatter to K_ii and K_b
                for i = 1:4
                    iGlobal = nodeList[i]
                    iFree = d_globalToFree[iGlobal]
                    iBoundary = d_globalToBoundary[iGlobal]

                    # Scatter to K_ii (only for interior nodes)
                    if iFree != -1
                        for j = 1:4
                            jGlobal = nodeList[j]
                            jFree = d_globalToFree[jGlobal]
                            k_val = Ke_fine[(i-1)*4 + j]

                            if jFree != -1
                                # Interior-interior: add to K_ii
                                col_start = d_colptr[iFree]
                                col_end = d_colptr[iFree + 1] - 1

                                # Binary search for jFree in rowidx
                                k_pos = binary_search_sparse(d_rowidx, jFree, col_start, col_end)
                                if k_pos != -1
                                    CUDA.@atomic d_valK_ii[offset_ii + k_pos] += k_val
                                end
                            else
                                # Interior-boundary: add to btmp (RHS for static condensation)
                                # Offset for this block's btmp
                                offset_btmp = (iel - 1) * nfree

                                for ir = 1:numVectorsToSolve
                                    phi_val = d_phi[offset_nodes + jGlobal, ir]
                                    btmp_idx = offset_btmp + iFree
                                    CUDA.@atomic d_btmp[btmp_idx, ir] -= k_val * phi_val
                                end
                            end
                        end
                    end

                    # Scatter to K_b (only for boundary nodes)
                    if iBoundary != -1
                        for j = 1:4
                            jGlobal = nodeList[j]
                            jGlobal_offset = offset_nodes + jGlobal
                            k_val = Ke_fine[(i-1)*4 + j]

                            # Boundary-all: add to K_b
                            row_start = d_rowptr_b[iBoundary]
                            row_end = d_rowptr_b[iBoundary + 1] - 1

                            # Binary search for jGlobal_offset in colidx
                            k_pos = binary_search_sparse(d_colidx_b, jGlobal_offset, row_start, row_end)
                            if k_pos != -1
                                CUDA.@atomic d_valK_b[offset_b + k_pos] += k_val
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# CUDA kernel to build block diagonal sparsity pattern (optimized for full thread utilization)
function build_block_diagonal_sparsity_kernel!(d_colptr, d_rowidx,
                                                d_colptrLocal, d_rowidxLocal,
                                                nb, nfree, nnz_ii)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    # Each thread handles one colptr entry
    if idx < nb * (nfree + 1)
        k = idx ÷ (nfree + 1)  # Block index
        j = idx % (nfree + 1)   # Local column index
        offset_col = k * nfree
        offset_nnz = k * nnz_ii
        @inbounds d_colptr[offset_col + j] = d_colptrLocal[j + 1] + offset_nnz
    end

    # Each thread handles one rowidx entry (offset to avoid overlap)
    idx_row = idx - nb * (nfree + 1)
    if idx_row >= 0 && idx_row < nb * nnz_ii
        k = idx_row ÷ nnz_ii
        i = idx_row % nnz_ii
        offset_nnz = k * nnz_ii
        @inbounds d_rowidx[offset_nnz + i + 1] = d_rowidxLocal[i + 1] + k * nfree
    end

    return nothing
end

# CUDA kernel to build block diagonal CSR sparsity pattern for K_b (optimized for full thread utilization)
function build_block_diagonal_csr_sparsity_kernel!(d_rowptr, d_colidx,
                                                    d_rowptrLocal, d_colidxLocal,
                                                    nb, nboundary, nnz_b, numNodes)
    # Global thread index
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    # Each thread handles one rowptr entry
    if idx < nb * (nboundary + 1)
        k = idx ÷ (nboundary + 1)  # Block index
        i = idx % (nboundary + 1)   # Local row index
        offset_row = k * nboundary
        offset_nnz = k * nnz_b
        @inbounds d_rowptr[offset_row + i] = d_rowptrLocal[i + 1] + offset_nnz
    end

    # Each thread handles one colidx entry (offset to avoid overlap)
    idx_col = idx - nb * (nboundary + 1)
    if idx_col >= 0 && idx_col < nb * nnz_b
        k = idx_col ÷ nnz_b
        i = idx_col % nnz_b
        offset_nnz = k * nnz_b
        offset_col = k * numNodes
        @inbounds d_colidx[offset_nnz + i + 1] = d_colidxLocal[i + 1] + offset_col
    end

    return nothing
end

# CUDA kernel to initialize d_utmp from Q1 basis functions
function init_utmp_kernel!(d_utmp, d_phiLocal, d_freeToGlobal,
                           nb, nfree, numVectorsToSolve)
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
function transfer_utmp_to_phi_kernel!(d_phi, d_utmp, d_freeToGlobal,
                                      nb, nfree, numNodes, numVectorsToSolve)
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
            d_phi[phi_offset + i_global, ir] = d_utmp[utmp_offset + i_local, ir]
        end

        # Compute 4th basis function from partition of unity: phi_4 = 1 - phi_1 - phi_2 - phi_3
        d_phi[phi_offset + i_global, 4] = 1.0 - d_phi[phi_offset + i_global, 1] -
                                                 d_phi[phi_offset + i_global, 2] -
                                                 d_phi[phi_offset + i_global, 3]
    end
    return nothing
end

# CUDA kernel to reconstruct fine-scale solution from coarse solution and basis functions
function reconstruct_fine_solution_kernel!(d_fine_x, d_fine_y, d_fine_u, d_coarse_solution,
                                           d_phi, d_cell_to_node, nx_coarse, ny_coarse,
                                           ratio, numNodes)
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

# CUDA kernel to compute coarse element RHS and stiffness:
# fe_coarse[iel, ir] = phi[iel, :, ir]' * rhs_fine[iel, :]
# Ke_coarse[iel, ir, jr] = phi[iel, :, ir]' * K_b[iel, :, :] * phi[iel, :, jr]
# Uses shared memory reduction for efficiency and to minimize atomic contention
function compute_fe_Ke_coarse_kernel!(d_fe_coarse, d_Ke_coarse, d_phi, d_rhs_fine,
                                       d_valK_b, d_rowptr_b, d_colidx_b, d_boundaryToGlobal,
                                       nb, numNodes, numVectors, nboundary, nnz_b)
    # Each block handles one coarse element
    iel = blockIdx().x
    tid = threadIdx().x

    if iel > nb
        return
    end

    # Offsets for this element's data
    offset_nodes = (iel - 1) * numNodes
    offset_b = (iel - 1) * nnz_b

    # Shared memory layout:
    # - First blockDim().x * numVectors for fe_coarse reduction
    # - Next blockDim().x * 16 for Ke_coarse local accumulation (16 = 4x4 matrix)
    sdata_fe = @cuDynamicSharedMem(Float64, blockDim().x * numVectors)
    sdata_Ke = @cuDynamicSharedMem(Float64, blockDim().x * 16, blockDim().x * numVectors * sizeof(Float64))

    # Initialize Ke_coarse local accumulator
    for i = 1:16
        @inbounds sdata_Ke[(tid - 1) * 16 + i] = 0.0
    end

    # ========================================================================
    # Part 1: Compute fe_coarse using reduction
    # ========================================================================

    # Each thread computes partial sums for all 4 basis functions
    for ir = 1:numVectors
        local_sum = 0.0
        # Stride loop: each thread processes multiple nodes if numNodes > blockDim().x
        i = tid
        while i <= numNodes
            @inbounds local_sum += d_phi[offset_nodes + i, ir] * d_rhs_fine[offset_nodes + i]
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

    # ========================================================================
    # Part 2: Compute Ke_coarse = phi^T * K_b * phi (using local accumulation)
    # ========================================================================

    # Parallelize over boundary rows: each thread handles some boundary rows
    # Accumulate locally in shared memory to avoid excessive atomics

    for iBoundary_local = tid:blockDim().x:nboundary
        if iBoundary_local <= nboundary
            # Global node index for this boundary node
            @inbounds iGlobal = d_boundaryToGlobal[iBoundary_local]

            # Row range in K_b for this boundary node
            @inbounds row_start = d_rowptr_b[iBoundary_local]
            @inbounds row_end = d_rowptr_b[iBoundary_local + 1] - 1

            # Iterate over non-zeros in this row
            for k_pos = row_start:row_end
                @inbounds jGlobal_offset = d_colidx_b[offset_b + k_pos]
                @inbounds k_val = d_valK_b[offset_b + k_pos]

                # Convert global offset to local node index
                jGlobal = jGlobal_offset - offset_nodes

                # Accumulate contribution to Ke_coarse for all basis function pairs
                # Store in thread-local shared memory (NO atomics here!)
                for ir = 1:numVectors
                    @inbounds phi_i = d_phi[offset_nodes + iGlobal, ir]
                    phi_i_k = phi_i * k_val

                    for jr = 1:numVectors
                        @inbounds phi_j = d_phi[offset_nodes + jGlobal, jr]
                        contribution = phi_i_k * phi_j

                        # Linear index into local Ke_coarse (4x4 = 16 elements)
                        ke_idx = (ir - 1) * numVectors + jr
                        # Accumulate in shared memory (thread-private)
                        @inbounds sdata_Ke[(tid - 1) * 16 + ke_idx] += contribution
                    end
                end
            end
        end
    end

    sync_threads()

    # Final reduction: sum all thread-local Ke_coarse contributions
    # Only use atomics at the very end (16 atomics per thread instead of thousands!)
    for i = 1:16
        local_val = sdata_Ke[(tid - 1) * 16 + i]
        if local_val != 0.0
            CUDA.@atomic d_Ke_coarse[iel, i] += local_val
        end
    end

    return nothing
end

# ============================================================================
# Main Program
# ============================================================================

function main()
    println("="^70)
    println("Julia 2D MFEM Assembly (Cuda with Coloring)")
    println("="^70)
    println()

    # ========================================================================
    # Mesh generation
    # ========================================================================

    nx, ny = 16, 16  # Coarse mesh (MFEM)
    ratio = 32        # Fine grid refinement ratio per coarse element

    println("Generating coarse mesh: $nx x $ny MFEM elements")
    println("MFEM ratio: $ratio (each coarse element has $(ratio)x$(ratio) fine elements)")
    println("Effective fine resolution: $(nx*ratio) x $(ny*ratio)")
    println()

    t0 = time()
    mesh = generate_mesh(nx, ny)
    mesh_time = time() - t0

    println("  Number of coarse elements: ", size(mesh.cell_to_node, 1))
    println("  Number of coarse nodes:    ", length(mesh.vertex_x))
    println("  Mesh generation:           ", @sprintf("%.2f ms", mesh_time * 1000))
    println()

    # ========================================================================
    # Build connectivity and coloring
    # ========================================================================

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
    println("  Matrix size: ", length(mesh.vertex_x), " x ", length(mesh.vertex_x))
    println("  Non-zeros:   ", nnz)
    println("  Number of colors: ", length(e2e_colors))
    for ic = 1:length(e2e_colors)
        println("    Color $ic: ", length(e2e_colors[ic]), " elements")
    end
    println("  Coarse Connectivity time: ", @sprintf("%.2f ms", conn_time * 1000))
    println()

    # Build the connectivty for one coarse element

    println("Building local mesh connectivity...")
    t0 = time()
    mesh_local = generate_mesh(ratio, ratio)
    n2n_row_ptr_local, n2n_col_idx_local, _ = build_mesh_connectivity(mesh_local)
    local_conn_time = time() - t0

    nnz = length(n2n_col_idx_local)
    println("  Local Matrix size: ", length(mesh_local.vertex_x), " x ", length(mesh_local.vertex_x))
    println("  Non-zeros:   ", nnz)
#    println("  Number of colors: ", length(e2e_colors))
#    for ic = 1:length(e2e_colors)
#        println("    Color $ic: ", length(e2e_colors[ic]), " elements")
#    end
    println("  Local Connectivity time: ", @sprintf("%.2f ms", local_conn_time * 1000))
    println()

  workspace = MFEWorkspace{Float64, 2, 4}(ratio);

  # Build DOF mappings and sparsity patterns using reusable kernel from FEMBase
  nfree, nboundary, nnz_ii, nnz_b = build_mfem_local_sparsity!(workspace, ratio, ratio)

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
          xi = Float64(ix) / ratio
          eta = Float64(iy) / ratio

          # Evaluate Q1 basis functions at (xi, eta)
          # Q1 nodes: [1] bottom-left, [2] bottom-right, [3] top-right, [4] top-left
          @inbounds phi[nodeID, 1] = (1.0 - xi) * (1.0 - eta)  # Bottom-left
          @inbounds phi[nodeID, 2] = xi * (1.0 - eta)          # Bottom-right
          @inbounds phi[nodeID, 3] = xi * eta                  # Top-right
          @inbounds phi[nodeID, 4] = (1.0 - xi) * eta          # Top-left
      end
  end

  nb = nx * ny  # Number of blocks

  d_phiLocal = cu(phi);
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
  @cuda threads=threads_per_block blocks=num_blocks build_block_diagonal_sparsity_kernel!(
      d_colptr_ii, d_rowidx_ii, d_colptrLocal, d_rowidxLocal, nb, nfree, nnz_ii
  )
  CUDA.synchronize()

  d_valK_ii = CUDA.zeros(Float64, nb * nnz_ii)

  # Build block diagonal CSR sparsity pattern for K_b (boundary-all coupling)
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
  @cuda threads=threads_per_block blocks=num_blocks build_block_diagonal_csr_sparsity_kernel!(
      d_rowptr_b, d_colidx_b, d_rowptrLocal_b, d_colidxLocal_b, nb, nboundary, nnz_b, numNodes
  )
  CUDA.synchronize()

  d_valK_b = CUDA.zeros(Float64, nb * nnz_b)

  # Allocate btmp (RHS for static condensation): (nb*nfree) × numVectorsToSolve
  numVectorsToSolve = 3
  d_btmp = CUDA.zeros(Float64, nb * nfree, numVectorsToSolve)

  # Upload globalToFree and globalToBoundary mappings to device
  d_globalToFree = cu(workspace.globalToFree)
  d_globalToBoundary = cu(workspace.globalToBoundary)

  # Allocate fine-grid RHS for all nodes (interior + boundary) in all coarse elements
  numNodes = (ratio + 1) * (ratio + 1)
  d_rhs_fine = CUDA.zeros(Float64, nb * numNodes)

  # Launch kernel to assemble K_ii, K_b, btmp, and rhs_fine
  println("GPU fine-grid assembly...")
  t0 = time()
  threads_per_block = 256
  num_blocks = cld(nb, threads_per_block)
  @cuda threads=threads_per_block blocks=num_blocks assemble_Kii_kernel!(
      d_valK_ii, d_btmp, d_rhs_fine, d_phi, d_colptr_ii, d_rowidx_ii, d_globalToFree,
      d_valK_b, d_rowptr_b, d_colidx_b, d_globalToBoundary,
      nx, ny, ratio, nfree, nnz_ii, nboundary, nnz_b, numVectorsToSolve, numNodes
  )
  CUDA.synchronize()
  gpu_fine_assembly_time = time() - t0
  println("  GPU fine assembly time: ", @sprintf("%.2f ms", gpu_fine_assembly_time * 1000))
  println()

  # Create the solution vector and initial guess
  d_utmp = CUDA.zeros(Float64, nb * nfree, numVectorsToSolve)
  d_freeToGlobal = cu(workspace.freeToGlobal)

  # Build initial guess from Q1 basis functions (partition of unity) - on device
  threads_per_block = 256
  num_blocks = cld(nb * nfree, threads_per_block)
  @cuda threads=threads_per_block blocks=num_blocks init_utmp_kernel!(
      d_utmp, d_phiLocal, d_freeToGlobal, nb, nfree, numVectorsToSolve
  )
  CUDA.synchronize()

  # Build the block diagonal sparse matrix K_ii on device
  n_total = nb * nfree
  K_ii_gpu = CuSparseMatrixCSC(d_colptr_ii, d_rowidx_ii, d_valK_ii, (n_total, n_total))

  # Solve K_ii * d_utmp = d_btmp using block CG
  println("Solving block CG on GPU...")
  t0 = time()
  # d_utmp contains initial guess and will contain solution after solve
  d_utmp, stats = block_cg(K_ii_gpu, d_btmp, d_utmp;
                           atol=1e-24, rtol=1e-12, itmax=1000, verbose=0)
  gpu_cg_time = time() - t0

  println("  Block CG converged: ", stats.solved)
  println("  Iterations: ", stats.niter)
  println("  Residual norm: ", stats.residuals[end])
  println("  GPU CG time: ", @sprintf("%.2f ms", gpu_cg_time * 1000))
  println()

  # Transfer solution from d_utmp back to d_phi
  threads_per_block = 256
  num_blocks = cld(nb * nfree, threads_per_block)
  @cuda threads=threads_per_block blocks=num_blocks transfer_utmp_to_phi_kernel!(
      d_phi, d_utmp, d_freeToGlobal, nb, nfree, numNodes, numVectorsToSolve
  )
  CUDA.synchronize()

  # Compute coarse element RHS and stiffness matrices
  # fe_coarse[iel, ir] = phi[iel, :, ir]' * rhs_fine[iel, :]
  # Ke_coarse[iel, ir*4+jr] = phi[iel, :, ir]' * K_b[iel, :, :] * phi[iel, :, jr]
  println("GPU coarse element assembly...")
  t0 = time()
  numVectors = 4
  d_fe_coarse = CUDA.zeros(Float64, nb, numVectors)
  d_Ke_coarse = CUDA.zeros(Float64, nb, numVectors * numVectors)

  # Upload boundaryToGlobal mapping to device
  d_boundaryToGlobal = cu(workspace.boundaryToGlobal)

  threads_per_block = 128  # Use 128 threads for better occupancy with shared memory
  num_blocks = nb  # One block per coarse element
  # Shared memory: fe_coarse reduction (threads * 4) + Ke_coarse local accumulation (threads * 16)
  shared_mem_size = threads_per_block * (numVectors + 16) * sizeof(Float64)

  @cuda threads=threads_per_block blocks=num_blocks shmem=shared_mem_size compute_fe_Ke_coarse_kernel!(
      d_fe_coarse, d_Ke_coarse, d_phi, d_rhs_fine,
      d_valK_b, d_rowptr_b, d_colidx_b, d_boundaryToGlobal,
      nb, numNodes, numVectors, nboundary, nnz_b
  )
  CUDA.synchronize()
  gpu_coarse_assembly_time = time() - t0
  println("  GPU coarse assembly time: ", @sprintf("%.2f ms", gpu_coarse_assembly_time * 1000))
  println()

# Create host copies of your device quantities
h_fe_coarse = Array(d_fe_coarse)
h_Ke_coarse = Array(d_Ke_coarse)

    # ========================================================================
    # Assembly for the coarse problem
    # ========================================================================

    println("="^70)
    println("Starting MFEM Assembly")
    println("="^70)

mat_values = zeros(Float64, length(n2n_col_idx))
rhs = zeros(Float64, length(mesh.vertex_x))

# Assemble global matrix and RHS from coarse element matrices using coloring
println("Assembling coarse system from element matrices...")
t0 = time()

# Loop over colors to assemble without race conditions
for (ic, elements) in enumerate(e2e_colors)
    # All elements in this color can be assembled in parallel (no conflicts)
    Threads.@threads for iel in elements
        # Get element nodes (Q1 element has 4 nodes)
        nodes = mesh.cell_to_node[iel, :]

        # Scatter RHS: rhs[gi] += fe_coarse[iel, i]
        for i = 1:4
            @inbounds gi = nodes[i]
            @inbounds rhs[gi] += h_fe_coarse[iel, i]
        end

        # Scatter stiffness: mat_values[k] += Ke_coarse[iel, i, j]
        for i = 1:4
            @inbounds gi = nodes[i]
            for j = 1:4
                @inbounds gj = nodes[j]
                # Find position in sparse matrix (uses binary search)
                k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                # h_Ke_coarse is stored as nb x 16 (flattened 4x4 matrices)
                ke_idx = (i - 1) * 4 + j
                @inbounds mat_values[k] += h_Ke_coarse[iel, ke_idx]
            end
        end
    end
end

assembly_time = time() - t0
println("  Assembly time: ", @sprintf("%.2f ms", assembly_time * 1000))
println()

    # ========================================================================
    # Apply boundary conditions
    # ========================================================================

    println("Applying boundary conditions...")
    t0 = time()
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global =
        apply_boundary_conditions(n2n_row_ptr, n2n_col_idx, mat_values, rhs, mesh.boundary_nodes)
    bc_time = time() - t0

    nfree = length(free_to_global)
    println("  Total DOFs:    ", length(mesh.vertex_x))
    println("  Boundary DOFs: ", length(mesh.boundary_nodes))
    println("  Free DOFs:     ", nfree)
    println("  Reduced nnz:   ", length(reduced_col_idx))
    println("  BC time:       ", @sprintf("%.2f ms", bc_time * 1000))
    println()

    # ========================================================================
    # Solve
    # ========================================================================

    println("Solving linear system...")

    # Build sparse matrix
    K = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)

    t0 = time()
    sol_free = K \ reduced_rhs
    solve_time = time() - t0

    println("  Solve time: ", @sprintf("%.2f ms", solve_time * 1000))
    println()

    # ========================================================================
    # Reconstruct Fine-Scale Solution (on GPU)
    # ========================================================================

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
    @cuda threads=threads_per_block blocks=num_blocks reconstruct_fine_solution_kernel!(
        d_fine_x, d_fine_y, d_fine_u, d_coarse_solution,
        d_phi, d_cell_to_node, nx, ny, ratio, numNodes
    )
    CUDA.synchronize()
    reconstruct_time = time() - t0

    # Transfer fine solution from GPU to CPU (only final result)
    fine_x = Array(d_fine_x)
    fine_y = Array(d_fine_y)
    fine_u = Array(d_fine_u)

    println("  Fine mesh size: ", length(fine_u), " nodes")
    println("  GPU reconstruction time: ", @sprintf("%.2f ms", reconstruct_time * 1000))
    println()

    # ========================================================================
    # Output
    # ========================================================================

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

    # ========================================================================
    # Performance summary
    # ========================================================================

    println("="^70)
    println("Performance Summary")
    println("="^70)
    println("CPU Timings:")
    println("  Mesh generation:     ", @sprintf("%8.2f ms", mesh_time * 1000))
    println("  Coarse connectivity: ", @sprintf("%8.2f ms", conn_time * 1000))
    println("  Local connectivity:  ", @sprintf("%8.2f ms", local_conn_time * 1000))
    println("  Coarse scatter:      ", @sprintf("%8.2f ms", assembly_time * 1000))
    println("  Boundary conditions: ", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Coarse solve:        ", @sprintf("%8.2f ms", solve_time * 1000))
    println()
    println("GPU Timings:")
    println("  Fine assembly:       ", @sprintf("%8.2f ms", gpu_fine_assembly_time * 1000))
    println("  Block CG solve:      ", @sprintf("%8.2f ms", gpu_cg_time * 1000))
    println("  Coarse assembly:     ", @sprintf("%8.2f ms", gpu_coarse_assembly_time * 1000))
    println("  Reconstruction:      ", @sprintf("%8.2f ms", reconstruct_time * 1000))
    println()
    println("  " * "-"^68)
    cpu_time = mesh_time + conn_time + local_conn_time + assembly_time + bc_time + solve_time
    gpu_time = gpu_fine_assembly_time + gpu_cg_time + gpu_coarse_assembly_time + reconstruct_time
    total_time = cpu_time + gpu_time
    println("  CPU Total:           ", @sprintf("%8.2f ms", cpu_time * 1000))
    println("  GPU Total:           ", @sprintf("%8.2f ms", gpu_time * 1000))
    println("  Total:               ", @sprintf("%8.2f ms", total_time * 1000))
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
