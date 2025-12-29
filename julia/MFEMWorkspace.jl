#!/usr/bin/env julia
#
# @file MFEMWorkspace.jl
# @brief Workspace structure for MFEM (Multiscale Finite Element Method) assembly
#
# This module provides pre-allocated workspace for MFEM element assembly
# to avoid allocations during the assembly loop.
#

module MFEMWorkspace

using LinearAlgebra

# Use FEMBase and PCG modules that are already loaded at the Main level
# (they must be included before this module is included)
using Main.FEMBase
using Main.PCG

export MFEWorkspace

"""
Workspace for MFEM element assembly (reused across elements to avoid allocations)

This workspace contains all pre-allocated arrays needed for MFEM static condensation:
- DOF mappings (interior vs boundary nodes)
- Sparse matrix arrays for K_ii (interior-interior) and K_b (boundary-all)
- RHS and basis function arrays
- PCG solver workspace for solving interior DOFs
- Element assembly workspace for fine-grid elements
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

    # Coarse element matrices (computed on GPU)
    Ke_coarse::Matrix{T}
    fe_coarse::Vector{T}

    # Coarse element data (separate from fine element workspace to avoid conflicts)
    nodes_coarse::Vector{Int}
    x_coarse::Vector{T}
    y_coarse::Vector{T}

    # PCG solver workspace (pre-allocated for in-place solving)
    pcg_workspace::PCG.PCGWorkspace{T}

    # Fine element assembly workspace (contains Ke, fe, x_coords, y_coords, nodes)
    element::FEMBase.ElementWorkspace{T, Dim, NNodes}
end

"""
    MFEWorkspace{T, Dim, NNodes}(ratio::Int) -> MFEWorkspace

Constructor for MFEWorkspace - pre-allocates all arrays based on refinement ratio.

# Arguments
- `ratio::Int`: Refinement ratio (fine grid has (ratio+1) × (ratio+1) nodes)

# Type Parameters
- `T`: Floating point type (typically Float64)
- `Dim`: Spatial dimension (typically 2)
- `NNodes`: Number of nodes per element (typically 4 for Q1)

# Returns
- Pre-allocated workspace with all arrays sized for the given ratio

# Notes
- Assumes square fine grid (same ratio in x and y)
- Pre-allocates sparse arrays with maximum possible sizes (9-point stencil)
- For ratio=32, creates workspace for 1089 fine nodes (33×33 grid)
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
        FEMBase.ElementWorkspace{T, Dim, NNodes}()
    )
end

end # module MFEMWorkspace
