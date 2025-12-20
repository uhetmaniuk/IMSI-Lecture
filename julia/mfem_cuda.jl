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

# ============================================================================
# MFEM-Specific Functions
# ============================================================================

"""
Workspace for MFEM element assembly (reused across elements to avoid allocations)
"""
struct MFEWorkspace{T<:AbstractFloat, Dim, NNodes}
    # DOF mapping arrays
    globalToFree::Vector{Int}
    freeToGlobal::Vector{Int}
    globalToBoundary::Vector{Int}
    boundaryToGlobal::Vector{Int}

    # K_ii (interior-interior) sparse matrix arrays
    colptr_ii::Vector{Int}
    rowidx_ii::Vector{Int}
    nzval_ii::Vector{T}

    # K_b (boundary-all) sparse matrix arrays
    matRowPtr_b::Vector{Int}
    matColIdx_b::Vector{Int}
    matValues_b::Vector{T}

    # RHS and basis function arrays
    rhs_fine::Vector{T}
    btmp::Matrix{T}
    phi::Matrix{T}
    utmp_init::Matrix{T}
    utmp::Matrix{T}

    # Coarse element matrices (output of compute_mfem_element!)
    Ke_coarse::Matrix{T}
    fe_coarse::Vector{T}

    # Coarse element data (separate from fine element workspace to avoid conflicts)
    nodes_coarse::Vector{Int}
    x_coarse::Vector{T}
    y_coarse::Vector{T}

    # PCG solver workspace (pre-allocated for in-place solving)
    pcg_workspace::PCG.PCGWorkspace{T}

    # Fine element assembly workspace (contains Ke, fe, x_coords, y_coords, nodes)
    element::ElementWorkspace{T, Dim, NNodes}
end

"""
Constructor for MFEWorkspace - pre-allocates all arrays based on ratio
"""
function MFEWorkspace{T, Dim, NNodes}(ratio::Int) where {T<:AbstractFloat, Dim, NNodes}
    numNodes = (ratio + 1) * (ratio + 1)
    numVectors = 4
    numVectorsToSolve = 3
    nfree = numNodes - 4 * ratio
    nboundary = 4 * ratio

    # Pre-allocate with maximum possible sizes (9-point stencil)
    max_nnz_ii = 9 * nfree
    max_nnz_b = 9 * nboundary

    return MFEWorkspace{T, Dim, NNodes}(
        # DOF mappings
        fill(-1, numNodes),           # globalToFree
        zeros(Int, nfree),            # freeToGlobal
        fill(-1, numNodes),           # globalToBoundary
        zeros(Int, nboundary),        # boundaryToGlobal

        # K_ii arrays
        zeros(Int, nfree + 1),        # colptr_ii
        zeros(Int, max_nnz_ii),       # rowidx_ii
        zeros(T, max_nnz_ii),         # nzval_ii

        # K_b arrays
        zeros(Int, nboundary + 1),    # matRowPtr_b
        zeros(Int, max_nnz_b),        # matColIdx_b
        zeros(T, max_nnz_b),          # matValues_b

        # RHS and basis functions
        zeros(T, numNodes),           # rhs_fine
        zeros(T, nfree, numVectorsToSolve),  # btmp
        zeros(T, numNodes, numVectors),      # phi
        zeros(T, nfree, numVectorsToSolve),  # utmp_init
        zeros(T, nfree, numVectorsToSolve),  # utmp

        # Coarse element matrices
        zeros(T, 4, 4),               # Ke_coarse
        zeros(T, 4),                  # fe_coarse

        # Coarse element data
        zeros(Int, 4),                # nodes_coarse
        zeros(T, 4),                  # x_coarse
        zeros(T, 4),                  # y_coarse

        # PCG solver workspace (pre-allocated for nfree interior DOFs)
        PCG.PCGWorkspace{T}(nfree),

        # Fine element workspace (contains Ke, fe, x_coords, y_coords, nodes)
        ElementWorkspace{T, Dim, NNodes}()
    )
end

"""
Compute MFEM basis functions and coarse element matrices
Returns: phi (basis functions on fine grid)
Uses optimized assembly with sparse K_ii (interior-interior) and K_b (boundary-all) matrices
Follows the C++ implementation in src/ScaledLaplacian.h (MFEMAssemblyFunctor::operator())
"""
function compute_mfem_element!(Ke_coarse::Matrix{Float64}, fe_coarse::Vector{Float64},
                              elem::Q1Element{Dim, NNodes}, workspace::MFEWorkspace{Float64, Dim, NNodes},
                              ratio::Int,
                              x_corner::Vector{Float64}, y_corner::Vector{Float64},
                              ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    numNodes = (ratio + 1) * (ratio + 1)
    numVectors = 4  # Number of coarse basis functions
    numVectorsToSolve = 3  # Solve only 3, get 4th from partition of unity

    # Use pre-allocated DOF mapping arrays from workspace
    globalToFree = workspace.globalToFree
    freeToGlobal = workspace.freeToGlobal
    globalToBoundary = workspace.globalToBoundary
    boundaryToGlobal = workspace.boundaryToGlobal

    # Reset DOF mappings
    fill!(globalToFree, -1)
    fill!(globalToBoundary, -1)

    nfree = 0;
    nboundary = 0;
    for iy = 0:ratio
        for ix = 0:ratio
            nodeID = ix + iy * (ratio + 1) + 1
            # Interior nodes
            if ix > 0 && ix < ratio && iy > 0 && iy < ratio
                nfree = nfree + 1
                @inbounds freeToGlobal[nfree] = nodeID
                @inbounds globalToFree[nodeID] = nfree;
            else
                # Boundary nodes
                nboundary = nboundary + 1;
                @inbounds boundaryToGlobal[nboundary] = nodeID
                @inbounds globalToBoundary[nodeID] = nboundary
            end
        end
    end

    # Use pre-allocated K_ii arrays from workspace
    colptr_ii = workspace.colptr_ii
    rowidx_ii = workspace.rowidx_ii
    nzval_ii = workspace.nzval_ii

    # Reset K_ii arrays
    fill!(nzval_ii, 0.0)

    # Build K_ii sparsity pattern (interior-interior coupling)
    # Note: We build the sparsity row-wise (CSR-style) but store in CSC format.
    # Since K_ii is symmetric, building A^T is equivalent to building A.
    colptr_ii[1] = 1

    for i = 1:nfree
        @inbounds iGlobal = freeToGlobal[i]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)

        # Count interior neighbors in 9-point stencil
        count = 1  # Diagonal
        hasWest = (ix > 1)
        hasEast = (ix < ratio - 1)
        hasSouth = (iy > 1)
        hasNorth = (iy < ratio - 1)

        # South neighbors
        if hasSouth
            count += hasWest + 1 + hasEast
        end
        # West and East
        count += hasWest + hasEast
        # North neighbors
        if hasNorth
            count += hasWest + 1 + hasEast
        end

        @inbounds colptr_ii[i + 1] = colptr_ii[i] + count
    end

    @inbounds nnz_ii = colptr_ii[nfree + 1] - 1

    # Fill K_ii sparsity (building row-wise, stored as CSC due to symmetry)
    for i = 1:nfree
        @inbounds iGlobal = freeToGlobal[i]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)
        @inbounds offset = colptr_ii[i]

        # Add interior neighbors
        # South row
        if iy > 1
            if ix > 1
                jGlobal = iGlobal - 1 - (ratio + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
            jGlobal = iGlobal - (ratio + 1)
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
            if ix < ratio - 1
                jGlobal = iGlobal + 1 - (ratio + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
        end
        # West
        if ix > 1
            jGlobal = iGlobal - 1
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
        end
        # Diagonal
        @inbounds rowidx_ii[offset] = i
        offset += 1
        # East
        if ix < ratio - 1
            jGlobal = iGlobal + 1
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
        end
        # North row
        if iy < ratio - 1
            if ix > 1
                jGlobal = iGlobal - 1 + (ratio + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
            jGlobal = iGlobal + (ratio + 1)
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
            if ix < ratio - 1
                jGlobal = iGlobal + 1 + (ratio + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
        end
    end

    # Use pre-allocated K_b arrays from workspace
    matRowPtr_b = workspace.matRowPtr_b
    matColIdx_b = workspace.matColIdx_b
    matValues_b = workspace.matValues_b

    # Reset K_b arrays
    fill!(matValues_b, 0.0)

    # Build K_b sparsity pattern (boundary-all coupling)
    @inbounds matRowPtr_b[1] = 1

    for i = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[i]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)

        # Count all neighbors (boundary and interior) in global numbering
        count = 1  # Diagonal
        hasWest = (ix > 0)
        hasEast = (ix < ratio)

        count += hasWest + hasEast
        if iy > 0
            count += 1 + hasWest + hasEast
        end
        if iy < ratio
            count += 1 + hasWest + hasEast
        end

        @inbounds matRowPtr_b[i + 1] = matRowPtr_b[i] + count
    end

    @inbounds nnz_b = matRowPtr_b[nboundary + 1] - 1

    # Fill K_b column indices (in global node numbering)
    for i = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[i]
        ix = (iGlobal - 1) % (ratio + 1)
        iy = (iGlobal - 1) ÷ (ratio + 1)
        @inbounds offset = matRowPtr_b[i]

        # Add all neighbors
        # South row
        if iy > 0
            if ix > 0
                @inbounds matColIdx_b[offset] = iGlobal - 1 - (ratio + 1)
                offset += 1
            end
            @inbounds matColIdx_b[offset] = iGlobal - (ratio + 1)
            offset += 1
            if ix < ratio
                @inbounds matColIdx_b[offset] = iGlobal + 1 - (ratio + 1)
                offset += 1
            end
        end
        # West
        if ix > 0
            @inbounds matColIdx_b[offset] = iGlobal - 1
            offset += 1
        end
        # Diagonal
        @inbounds matColIdx_b[offset] = iGlobal
        offset += 1
        # East
        if ix < ratio
            @inbounds matColIdx_b[offset] = iGlobal + 1
            offset += 1
        end
        # North row
        if iy < ratio
            if ix > 0
                @inbounds matColIdx_b[offset] = iGlobal - 1 + (ratio + 1)
                offset += 1
            end
            @inbounds matColIdx_b[offset] = iGlobal + (ratio + 1)
            offset += 1
            if ix < ratio
                @inbounds matColIdx_b[offset] = iGlobal + 1 + (ratio + 1)
                offset += 1
            end
        end
    end

    # Use pre-allocated RHS and basis function arrays from workspace
    rhs_fine = workspace.rhs_fine
    btmp = workspace.btmp
    phi = workspace.phi
    utmp_init = workspace.utmp_init

    # Reset arrays
    fill!(rhs_fine, 0.0)
    fill!(btmp, 0.0)

    # Set corner values
    @inbounds phi[1, 1] = 1.0                                    # Bottom-left corner -> Basis 0
    @inbounds phi[ratio + 1, 2] = 1.0                           # Bottom-right corner -> Basis 1
    @inbounds phi[(ratio + 1) * (ratio + 1), 3] = 1.0           # Top-right corner -> Basis 2
    @inbounds phi[(ratio + 1) * ratio + 1, 4] = 1.0             # Top-left corner -> Basis 3

    # Set edge values (Q1 interpolation), skip corners
    for is = 1:ratio-1
        s = Float64(is) / ratio
        # Left edge (ix = 0)
        nodeID = is * (ratio + 1) + 1
        @inbounds phi[nodeID, 1] = 1.0 - s  # Basis 0
        @inbounds phi[nodeID, 4] = s        # Basis 3
        # Right edge (ix = ratio)
        nodeID = ratio + is * (ratio + 1) + 1
        @inbounds phi[nodeID, 2] = 1.0 - s  # Basis 1
        @inbounds phi[nodeID, 3] = s        # Basis 2
        # Bottom edge (iy = 0)
        nodeID = is + 1
        @inbounds phi[nodeID, 1] = 1.0 - s  # Basis 0
        @inbounds phi[nodeID, 2] = s        # Basis 1
        # Top edge (iy = ratio)
        nodeID = is + ratio * (ratio + 1) + 1
        @inbounds phi[nodeID, 4] = 1.0 - s  # Basis 3
        @inbounds phi[nodeID, 3] = s        # Basis 2
    end

    # Assemble fine elements with scatter to K_ii, K_b, and btmp
    @inbounds hx = (x_corner[2] - x_corner[1]) / ratio
    @inbounds hy = (y_corner[4] - y_corner[1]) / ratio

    # Use pre-allocated arrays from element workspace (no redundant allocation!)
    x_fine = workspace.element.x_coords
    y_fine = workspace.element.y_coords
    nodeList = workspace.element.nodes

    for iy = 0:ratio-1
        for ix = 0:ratio-1
            # Fine element corners
            @inbounds x_fine[1] = x_corner[1] + ix * hx
            @inbounds y_fine[1] = y_corner[1] + iy * hy
            @inbounds x_fine[2] = x_fine[1] + hx
            @inbounds y_fine[2] = y_fine[1]
            @inbounds x_fine[3] = x_fine[2]
            @inbounds y_fine[3] = y_fine[2] + hy
            @inbounds x_fine[4] = x_fine[1]
            @inbounds y_fine[4] = y_fine[3]

            # Assemble fine element (Ke and fe are in workspace.element)
            assemble_element!(elem, workspace.element, x_fine, y_fine, ax_func, ay_func, f_func)

            # Extract Ke and fe from element workspace
            Ke_fine = workspace.element.Ke
            fe_fine = workspace.element.fe

            # Fine element node list
            nodeList[1] = ix + iy * (ratio + 1) + 1
            nodeList[2] = ix + 1 + iy * (ratio + 1) + 1
            nodeList[3] = ix + 1 + (iy + 1) * (ratio + 1) + 1
            nodeList[4] = ix + (iy + 1) * (ratio + 1) + 1

            # Scatter RHS
            for i = 1:4
                @inbounds gi = nodeList[i]
                @inbounds rhs_fine[gi] += fe_fine[i]
            end

            # Scatter stiffness to K_ii, K_b, and btmp (hoisted branches for performance)
            for i = 1:4
                @inbounds iGlobal = nodeList[i]
                @inbounds iFree = globalToFree[iGlobal]

                if iFree != -1
                    # Interior row - process all j together
                    for j = 1:4
                        @inbounds jGlobal = nodeList[j]
                        @inbounds jFree = globalToFree[jGlobal]
                        @inbounds k_val = Ke_fine[i, j]

                        if jFree != -1
                            # Interior-interior: add to K_ii
                            k = find_matrix_position(colptr_ii, rowidx_ii, iFree, jFree)
                            @inbounds nzval_ii[k] += k_val
                        else
                            # Interior-boundary: add to btmp (RHS for static condensation)
                            @simd for ir = 1:numVectorsToSolve
                                @inbounds phi_val = phi[jGlobal, ir]
                                @inbounds btmp[iFree, ir] -= k_val * phi_val
                            end
                        end
                    end
                else
                    # Boundary row - process all j together
                    @inbounds iBoundary = globalToBoundary[iGlobal]
                    for j = 1:4
                        @inbounds jGlobal = nodeList[j]
                        @inbounds k_val = Ke_fine[i, j]
                        k = find_matrix_position(matRowPtr_b, matColIdx_b, iBoundary, jGlobal)
                        @inbounds matValues_b[k] += k_val
                    end
                end
            end
        end
    end

    # Build initial guess for PCG based on Q1 shape functions
    for i = 1:nfree
        @inbounds gi = freeToGlobal[i]

        # Compute normalized coordinates (xi, eta) of this fine node within coarse element
        # Fine nodes are arranged in a (ratio+1) x (ratio+1) grid
        iy_fine = (gi - 1) ÷ (ratio + 1)
        ix_fine = (gi - 1) % (ratio + 1)
        xi = Float64(ix_fine) / ratio   # Normalized coordinate [0, 1]
        eta = Float64(iy_fine) / ratio  # Normalized coordinate [0, 1]

        # Evaluate Q1 shape functions at this position
        # Q1 nodes: [1] bottom-left, [2] bottom-right, [3] top-right, [4] top-left
        N1 = (1.0 - xi) * (1.0 - eta)  # Basis 1 (bottom-left)
        N2 = xi * (1.0 - eta)          # Basis 2 (bottom-right)
        N3 = xi * eta                  # Basis 3 (top-right)

        # Initial guess for interior values of first 3 basis functions
        @inbounds utmp_init[i, 1] = N1
        @inbounds utmp_init[i, 2] = N2
        @inbounds utmp_init[i, 3] = N3
    end

    # Solve for interior basis functions using in-place PCG with SSOR preconditioning
    # Use only the filled portions of the pre-allocated arrays
    K_ii = SparseMatrixCSC(nfree, nfree, colptr_ii[1:nfree+1], rowidx_ii[1:nnz_ii], nzval_ii[1:nnz_ii])

    # Copy initial guess to utmp (will be modified in-place)
    utmp = workspace.utmp
    @. utmp = utmp_init

    # Solve in-place (ZERO allocations!)
    pcg_info = pcg_solve_ssor!(utmp, workspace.pcg_workspace, K_ii, btmp,
                               omega=1.0, num_ssor_sweeps=1,
                               tol=1e-12, maxiter=1000, verbose=false)

    # Check convergence
    #if !pcg_info.converged
    #    @warn "PCG-SSOR did not converge: iterations=$(pcg_info.iterations), residual=$(pcg_info.residual_norm)"
    #end

    # Reconstruct full basis functions: copy interior solutions and compute 4th basis
    for i = 1:nfree
        @inbounds gi = freeToGlobal[i]
        sum_val = 0.0
        for ir = 1:numVectorsToSolve
            @inbounds phi[gi, ir] = utmp[i, ir]
            sum_val += utmp[i, ir]
        end
        # 4th basis function from partition of unity
        @inbounds phi[gi, 4] = 1.0 - sum_val
    end

    # Compute coarse element RHS: fe_coarse[ir] = phi[:,ir]' * rhs_fine
    for ir = 1:numVectors
        sum_val = 0.0
        @simd for i = 1:numNodes
            @inbounds sum_val += phi[i, ir] * rhs_fine[i]
        end
        @inbounds fe_coarse[ir] = sum_val
    end

    # Compute coarse element stiffness using Schur complement: Ke_coarse = phi^T * K * phi
    # This is computed efficiently using K_b (boundary rows only):
    # kele = phi_b^T * K_b * phi
    # where K_b includes boundary-interior and boundary-boundary coupling
    fill!(Ke_coarse, 0.0)

    for iBoundary = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[iBoundary]
        @inbounds rowBegin = matRowPtr_b[iBoundary]
        @inbounds rowEnd = matRowPtr_b[iBoundary + 1] - 1

        for k = rowBegin:rowEnd
            @inbounds jGlobal = matColIdx_b[k]  # K_b uses global node numbering for columns
            @inbounds k_val = matValues_b[k]

            # Accumulate contribution to Ke_coarse for all basis function pairs
            # Hoist phi_i * k_val out of jr loop for better performance
            for ir = 1:numVectors
                @inbounds phi_i_k = phi[iGlobal, ir] * k_val
                @simd for jr = 1:numVectors
                    @inbounds phi_j = phi[jGlobal, jr]
                    @inbounds Ke_coarse[ir, jr] += phi_i_k * phi_j
                end
            end
        end
    end

    return phi
end

# ============================================================================
# System Assembly with MFEM
# ============================================================================

"""
Assemble global system with MFEM elements (thread-parallel with coloring)
Returns: mat_values, rhs, basis_functions, compute_mfem_time
Uses optimized kernels from FEMBase with thread-local workspaces
"""
function assemble_system_mfem(mesh::Mesh, elem::Q1Element{Dim, NNodes}, ratio::Int,
                             n2n_row_ptr::Vector{Int}, n2n_col_idx::Vector{Int},
                             e2e_colors::Vector{Vector{Int}},
                             ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    nnodes = length(mesh.vertex_x)
    nnz = length(n2n_col_idx)
    nel = size(mesh.cell_to_node, 1)

    # Initialize global arrays
    mat_values = zeros(nnz)
    rhs = zeros(nnodes)

    # Store basis functions for each element (for fine-scale reconstruction)
    numFineNodes = (ratio + 1) * (ratio + 1)
    basis_functions = [zeros(numFineNodes, 4) for _ = 1:nel]

    # Thread-local workspaces (Ke_coarse and fe_coarse are now in workspace)
    max_tid = Threads.maxthreadid()
    workspace_local = [MFEWorkspace{Float64, Dim, NNodes}(ratio) for _ = 1:max_tid]

    # Loop over colors
    for (ic, elements) in enumerate(e2e_colors)
        # All elements in this color can be assembled in parallel (no conflicts)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            @inbounds workspace = workspace_local[tid]

            # Use pre-allocated coarse element arrays from workspace (avoid allocation!)
            nodes = workspace.nodes_coarse
            x = workspace.x_coarse
            y = workspace.y_coarse
            Ke = workspace.Ke_coarse
            fe = workspace.fe_coarse

            # Extract element data (manual copy, no allocation)
            @inbounds for i = 1:4
                nodes[i] = mesh.cell_to_node[iel, i]
            end
            @inbounds for i = 1:4
                ni = nodes[i]
                x[i] = mesh.vertex_x[ni]
                y[i] = mesh.vertex_y[ni]
            end

            # Assemble MFEM element and get basis functions using optimized kernel
            @inbounds basis_functions[iel] = compute_mfem_element!(Ke, fe, elem, workspace, ratio, x, y, ax_func, ay_func, f_func)

            # Scatter to global arrays (no race condition within same color)
            for i = 1:4
                @inbounds gi = nodes[i]
                @inbounds rhs[gi] += fe[i]
                for j = 1:4
                    @inbounds gj = nodes[j]
                    # Use binary search from FEMBase (assumes sorted columns)
                    k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                    @inbounds mat_values[k] += Ke[i, j]
                end
            end
        end
    end

    return mat_values, rhs, basis_functions
end

# ============================================================================
# Fine-Scale Solution Reconstruction
# ============================================================================

"""
Reconstruct fine-scale solution from coarse solution using MFEM basis functions
Uses graph coloring to avoid race conditions on shared fine nodes
"""
function reconstruct_fine_solution(mesh::Mesh, coarse_solution::Vector{Float64},
                                  basis_functions::Vector{Matrix{Float64}},
                                  e2e_colors::Vector{Vector{Int}}, ratio::Int)
    nel = size(mesh.cell_to_node, 1)
    nx_coarse = Int(sqrt(nel))  # Assuming square mesh
    ny_coarse = nx_coarse

    # Fine mesh dimensions
    nx_fine = nx_coarse * ratio
    ny_fine = ny_coarse * ratio
    nfine = (nx_fine + 1) * (ny_fine + 1)

    # Fine mesh coordinates and solution
    fine_x = zeros(nfine)
    fine_y = zeros(nfine)
    fine_u = zeros(nfine)

    # Build fine mesh (parallel - independent writes)
    Threads.@threads for j = 0:ny_fine
        for i = 0:nx_fine
            idx = (i + 1) + j * (nx_fine + 1)
            @inbounds fine_x[idx] = Float64(i) / nx_fine
            @inbounds fine_y[idx] = Float64(j) / ny_fine
        end
    end

    # Thread-local temporary arrays (avoid allocations in hot loop)
    max_tid = Threads.maxthreadid()
    nodes_local = [zeros(Int, 4) for _ = 1:max_tid]
    u_local = [zeros(4) for _ = 1:max_tid]

    # Reconstruct solution element by element (parallel with coloring to avoid races)
    for (ic, elements) in enumerate(e2e_colors)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            nodes = nodes_local[tid]
            u_coarse_elem = u_local[tid]

            # Coarse element indices
            iel_1d = iel - 1
            ix_coarse = iel_1d % nx_coarse
            iy_coarse = iel_1d ÷ nx_coarse

            # Extract coarse element nodes and solution values (no allocation!)
            @inbounds for i = 1:4
                nodes[i] = mesh.cell_to_node[iel, i]
                u_coarse_elem[i] = coarse_solution[nodes[i]]
            end

            # Basis functions for this element
            @inbounds phi = basis_functions[iel]

            # Reconstruct fine solution: u_fine = phi * u_coarse
            for iy_local = 0:ratio
                for ix_local = 0:ratio
                    # Fine node global indices
                    ix_fine = ix_coarse * ratio + ix_local
                    iy_fine = iy_coarse * ratio + iy_local
                    fine_node = ix_fine + iy_fine * (nx_fine + 1) + 1

                    # Local fine node index
                    local_node = ix_local + iy_local * (ratio + 1) + 1

                    # Interpolate: u_fine = sum_i phi[local_node, i] * u_coarse[i]
                    u_val = 0.0
                    @simd for i = 1:4
                        @inbounds u_val += phi[local_node, i] * u_coarse_elem[i]
                    end
                    @inbounds fine_u[fine_node] = u_val
                end
            end
        end
    end

    return fine_x, fine_y, fine_u
end

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

                                # Linear search for jFree in rowidx
                                for k_pos = col_start:col_end
                                    if d_rowidx[k_pos] == jFree
                                        CUDA.@atomic d_valK_ii[offset_ii + k_pos] += k_val
                                        break
                                    end
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

                            # Linear search for jGlobal_offset in colidx
                            for k_pos = row_start:row_end
                                if d_colidx_b[k_pos] == jGlobal_offset
                                    CUDA.@atomic d_valK_b[offset_b + k_pos] += k_val
                                    break
                                end
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
function build_block_diagonal_sparsity_kernel!(d_colptr, d_rowidx,
                                                d_colptrLocal, d_rowidxLocal,
                                                nb, nfree, nnz_ii)
    # One thread per block
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    if k < nb
        # Build colptr for block k
        offset_col = k * nfree
        offset_nnz = k * nnz_ii

        for j = 1:nfree+1
            d_colptr[offset_col + j] = d_colptrLocal[j] + offset_nnz
        end

        # Build rowidx for block k
        for i = 1:nnz_ii
            d_rowidx[offset_nnz + i] = d_rowidxLocal[i] + k * nfree
        end
    end
    return nothing
end

# CUDA kernel to build block diagonal CSR sparsity pattern for K_b
function build_block_diagonal_csr_sparsity_kernel!(d_rowptr, d_colidx,
                                                    d_rowptrLocal, d_colidxLocal,
                                                    nb, nboundary, nnz_b, numNodes)
    # One thread per block
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    if k < nb
        # Build rowptr for block k
        offset_row = k * nboundary
        offset_nnz = k * nnz_b

        for i = 1:nboundary+1
            d_rowptr[offset_row + i] = d_rowptrLocal[i] + offset_nnz
        end

        # Build colidx for block k (column indices are in global numbering within element)
        offset_col = k * numNodes
        for i = 1:nnz_b
            d_colidx[offset_nnz + i] = d_colidxLocal[i] + offset_col
        end
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

# CUDA kernel to compute coarse element RHS and stiffness:
# fe_coarse[iel, ir] = phi[iel, :, ir]' * rhs_fine[iel, :]
# Ke_coarse[iel, ir, jr] = phi[iel, :, ir]' * K_b[iel, :, :] * phi[iel, :, jr]
# Uses shared memory reduction for efficiency
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

    # Shared memory for partial sums (4 basis functions) for fe_coarse
    sdata = @cuDynamicSharedMem(Float64, (blockDim().x, numVectors))

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
        @inbounds sdata[tid, ir] = local_sum
    end

    sync_threads()

    # Reduction in shared memory
    s = blockDim().x ÷ 2
    while s > 0
        if tid <= s
            for ir = 1:numVectors
                @inbounds sdata[tid, ir] += sdata[tid + s, ir]
            end
        end
        sync_threads()
        s = s ÷ 2
    end

    # Thread 1 writes fe_coarse result
    if tid == 1
        for ir = 1:numVectors
            @inbounds d_fe_coarse[iel, ir] = sdata[1, ir]
        end
    end

    sync_threads()

    # ========================================================================
    # Part 2: Compute Ke_coarse = phi^T * K_b * phi
    # ========================================================================

    # Parallelize over boundary rows: each thread handles some boundary rows
    # Since Ke_coarse is small (4x4), we use atomic adds for simplicity

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
                for ir = 1:numVectors
                    @inbounds phi_i = d_phi[offset_nodes + iGlobal, ir]
                    phi_i_k = phi_i * k_val

                    for jr = 1:numVectors
                        @inbounds phi_j = d_phi[offset_nodes + jGlobal, jr]
                        contribution = phi_i_k * phi_j

                        # Linear index into Ke_coarse (stored as nb x 16 matrix)
                        ke_idx = (ir - 1) * numVectors + jr
                        CUDA.@atomic d_Ke_coarse[iel, ke_idx] += contribution
                    end
                end
            end
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
    # Problem setup
    # ========================================================================

    # Material coefficients and forcing term (varying coefficients)
    ax(x, y, z) = 1.0
    ay(x, y, z) = 1.0
    function f(x, y, z)
        # Manufactured solution: u = sin(π*x) * sin(π*y)
        return 2.0 * π^2 * sin(π * x) * sin(π * y)
    end

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
    n2n_row_ptr, n2n_col_idx, _ = build_mesh_connectivity(mesh_local)
    conn_time = time() - t0

    nnz = length(n2n_col_idx)
    println("  Local Matrix size: ", length(mesh_local.vertex_x), " x ", length(mesh_local.vertex_x))
    println("  Non-zeros:   ", nnz)
#    println("  Number of colors: ", length(e2e_colors))
#    for ic = 1:length(e2e_colors)
#        println("    Color $ic: ", length(e2e_colors[ic]), " elements")
#    end
    println("  Local Connectivity time: ", @sprintf("%.2f ms", conn_time * 1000))
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

  threads_per_block = 256
  num_blocks = cld(nb, threads_per_block)
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

  threads_per_block = 256
  num_blocks = cld(nb, threads_per_block)
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
  threads_per_block = 256
  num_blocks = cld(nb, threads_per_block)
  @cuda threads=threads_per_block blocks=num_blocks assemble_Kii_kernel!(
      d_valK_ii, d_btmp, d_rhs_fine, d_phi, d_colptr_ii, d_rowidx_ii, d_globalToFree,
      d_valK_b, d_rowptr_b, d_colidx_b, d_globalToBoundary,
      nx, ny, ratio, nfree, nnz_ii, nboundary, nnz_b, numVectorsToSolve, numNodes
  )

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
  # d_utmp contains initial guess and will contain solution after solve
  d_utmp, stats = block_cg(K_ii_gpu, d_btmp, d_utmp;
                           atol=1e-24, rtol=1e-12, itmax=1000, verbose=0)

  println("Block CG converged: ", stats.solved)
  println("  Iterations: ", stats.niter)
  println("  Residual norm: ", stats.residuals[end])

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
  numVectors = 4
  d_fe_coarse = CUDA.zeros(Float64, nb, numVectors)
  d_Ke_coarse = CUDA.zeros(Float64, nb, numVectors * numVectors)

  # Upload boundaryToGlobal mapping to device
  d_boundaryToGlobal = cu(workspace.boundaryToGlobal)

  threads_per_block = 128  # Use 128 threads for better occupancy with shared memory
  num_blocks = nb  # One block per coarse element
  shared_mem_size = threads_per_block * numVectors * sizeof(Float64)

  @cuda threads=threads_per_block blocks=num_blocks shmem=shared_mem_size compute_fe_Ke_coarse_kernel!(
      d_fe_coarse, d_Ke_coarse, d_phi, d_rhs_fine,
      d_valK_b, d_rowptr_b, d_colidx_b, d_boundaryToGlobal,
      nb, numNodes, numVectors, nboundary, nnz_b
  )
  CUDA.synchronize()

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
    # Reconstruct Fine-Scale Solution
    # ========================================================================

    # ========================================================================
    # Output
    # ========================================================================

    # ========================================================================
    # Performance summary
    # ========================================================================

end

# Run main program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
