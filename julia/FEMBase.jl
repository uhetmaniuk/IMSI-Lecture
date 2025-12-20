#!/usr/bin/env julia
#
# @file FEMBase.jl
# @brief Shared base kernels for FEM assembly in Julia
#
# This module contains optimized, reusable components for finite element
# assembly that can be shared between different FEM implementations.
#

module FEMBase

using SparseArrays
using LinearAlgebra
using Printf

export Q1Element, ElementWorkspace
export Mesh, generate_mesh
export transpose_graph, combine_graphs, color_graph
export build_mesh_connectivity
export assemble_element!
export find_matrix_position, apply_boundary_conditions
export write_solution
export shape_functions!, shape_gradients!
export build_mfem_local_sparsity!

# ============================================================================
# Q1 Finite Element (2D Bilinear Quadrilateral)
# ============================================================================

"""
Q1 shape functions and gradients on reference element [-1,1]²
Type parameters Dim and NNodes are compile-time constants for optimization
"""
struct Q1Element{Dim, NNodes}
    # Gauss quadrature points and weights for 2x2 rule
    qpts::Vector{Tuple{Float64, Float64}}
    qwts::Vector{Float64}

    function Q1Element{Dim, NNodes}() where {Dim, NNodes}
        # 2x2 Gauss quadrature
        gp = 1.0 / sqrt(3.0)
        qpts = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
        qwts = [1.0, 1.0, 1.0, 1.0]
        new{Dim, NNodes}(qpts, qwts)
    end
end

# Convenience constructor
Q1Element() = Q1Element{2, 4}()

# Accessor methods for dimension and node count
Base.getproperty(::Q1Element{Dim, NNodes}, ::Val{:dim}) where {Dim, NNodes} = Dim
Base.getproperty(::Q1Element{Dim, NNodes}, ::Val{:numNodes}) where {Dim, NNodes} = NNodes
function Base.getproperty(elem::Q1Element, name::Symbol)
    if name === :dim
        return getproperty(elem, Val(:dim))
    elseif name === :numNodes
        return getproperty(elem, Val(:numNodes))
    else
        return getfield(elem, name)
    end
end

"""
Evaluate Q1 shape functions at (ξ, η) - in-place version
Fills N with 4 shape function values [N1, N2, N3, N4]
Node ordering: (ξ,η) = [(-1,-1), (1,-1), (1,1), (-1,1)]
"""
function shape_functions!(N::Vector{Float64}, ξ::Float64, η::Float64)
    @inbounds N[1] = 0.25 * (1 - ξ) * (1 - η)
    @inbounds N[2] = 0.25 * (1 + ξ) * (1 - η)
    @inbounds N[3] = 0.25 * (1 + ξ) * (1 + η)
    @inbounds N[4] = 0.25 * (1 - ξ) * (1 + η)
    return nothing
end

"""
Evaluate Q1 shape function gradients in reference coordinates - in-place version
Fills dN (4x2 matrix) where row i contains [∂Ni/∂ξ, ∂Ni/∂η]
"""
function shape_gradients!(dN::Matrix{Float64}, ξ::Float64, η::Float64)
    # ∂N/∂ξ
    @inbounds dN[1, 1] = -0.25 * (1 - η)
    @inbounds dN[2, 1] =  0.25 * (1 - η)
    @inbounds dN[3, 1] =  0.25 * (1 + η)
    @inbounds dN[4, 1] = -0.25 * (1 + η)
    # ∂N/∂η
    @inbounds dN[1, 2] = -0.25 * (1 - ξ)
    @inbounds dN[2, 2] = -0.25 * (1 + ξ)
    @inbounds dN[3, 2] =  0.25 * (1 + ξ)
    @inbounds dN[4, 2] =  0.25 * (1 - ξ)
    return nothing
end

# ============================================================================
# Workspace for Element Assembly
# ============================================================================

"""
Workspace for element assembly (reused across elements to avoid allocations)
"""
struct ElementWorkspace{T<:AbstractFloat, Dim, NNodes}
    N::Vector{T}
    dN_ref::Matrix{T}
    dN_phys::Matrix{T}
    J::Matrix{T}
    invJ::Matrix{T}
    Ke::Matrix{T}
    fe::Vector{T}
    # Pre-allocated arrays for element data extraction
    nodes::Vector{Int}
    x_coords::Vector{T}
    y_coords::Vector{T}

    function ElementWorkspace{T, Dim, NNodes}() where {T<:AbstractFloat, Dim, NNodes}
        new{T, Dim, NNodes}(
            zeros(T, NNodes),
            zeros(T, NNodes, Dim),
            zeros(T, NNodes, Dim),
            zeros(T, Dim, Dim),
            zeros(T, Dim, Dim),
            zeros(T, NNodes, NNodes),
            zeros(T, NNodes),
            zeros(Int, NNodes),
            zeros(T, NNodes),
            zeros(T, NNodes)
        )
    end
end

# ============================================================================
# Mesh Generation
# ============================================================================

"""
Mesh data structure for 2D Q1 elements
"""
struct Mesh
    vertex_x::Vector{Float64}
    vertex_y::Vector{Float64}
    cell_to_node::Matrix{Int}  # nel x 4 (node indices, 1-based)
    boundary_nodes::Vector{Int}
end

"""
Generate structured 2D rectangular mesh with Q1 elements
Domain: [0,1] x [0,1]
Returns: Mesh object
"""
function generate_mesh(nx::Int, ny::Int)
    # Generate vertices
    nv = (nx + 1) * (ny + 1)
    vertex_x = zeros(nv)
    vertex_y = zeros(nv)

    idx = 1
    for j = 0:ny
        for i = 0:nx
            vertex_x[idx] = Float64(i) / nx
            vertex_y[idx] = Float64(j) / ny
            idx += 1
        end
    end

    # Generate elements (Q1 = 4 nodes per element)
    nel = nx * ny
    cell_to_node = zeros(Int, nel, 4)

    iel = 1
    for j = 0:ny-1
        for i = 0:nx-1
            n1 = j * (nx + 1) + i + 1      # Bottom-left
            n2 = n1 + 1                     # Bottom-right
            n3 = n2 + (nx + 1)              # Top-right
            n4 = n1 + (nx + 1)              # Top-left
            cell_to_node[iel, :] = [n1, n2, n3, n4]
            iel += 1
        end
    end

    # Identify boundary nodes
    boundary_nodes = Int[]
    for i = 1:nv
        x, y = vertex_x[i], vertex_y[i]
        if x ≈ 0.0 || x ≈ 1.0 || y ≈ 0.0 || y ≈ 1.0
            push!(boundary_nodes, i)
        end
    end

    return Mesh(vertex_x, vertex_y, cell_to_node, boundary_nodes)
end

# ============================================================================
# Graph Operations
# ============================================================================

"""
Transpose a graph stored as CSR (row_ptr, col_idx) - Parallel version
"""
function transpose_graph(n::Int, m::Int, row_ptr::Vector{Int}, col_idx::Vector{Int})
    nnz = length(col_idx)
    # Use maxthreadid() to get the maximum possible thread ID
    max_tid = Threads.maxthreadid()

    # Phase 1: Count entries per column (parallel reduction)
    col_counts_per_thread = [zeros(Int, m) for _ = 1:max_tid]

    Threads.@threads for i = 1:n
        tid = Threads.threadid()
        local_count = col_counts_per_thread[tid]
        for k = row_ptr[i]:row_ptr[i + 1] - 1
            local_count[col_idx[k]] += 1
        end
    end

    # Reduce counts from all threads
    col_count = zeros(Int, m)
    for tid = 1:max_tid
        col_count .+= col_counts_per_thread[tid]
    end

    # Phase 2: Build new row_ptr (parallel prefix sum)
    new_row_ptr = zeros(Int, m + 1)
    new_row_ptr[1] = 1
    for i = 1:m
        new_row_ptr[i + 1] = new_row_ptr[i] + col_count[i]
    end

    # Phase 3: Fill new_col_idx (parallel with atomics)
    new_col_idx = zeros(Int, nnz)
    col_offset = [Threads.Atomic{Int}(0) for _ = 1:m]

    # Parallel fill using atomic operations
    Threads.@threads for i = 1:n
        for k = row_ptr[i]:row_ptr[i + 1] - 1
            j = col_idx[k]
            # Atomically get and increment the offset for column j
            offset = Threads.atomic_add!(col_offset[j], 1)
            pos = new_row_ptr[j] + offset
            new_col_idx[pos] = i
        end
    end

    return new_row_ptr, new_col_idx
end

"""
Combine two graphs: A → B and B → C to get A → C - Parallel version
"""
function combine_graphs(na::Int, nb::Int, nc::Int,
                       ab_row_ptr::Vector{Int}, ab_col_idx::Vector{Int},
                       bc_row_ptr::Vector{Int}, bc_col_idx::Vector{Int})

    # Phase 1: Count unique entries
    ac_row_ptr = zeros(Int, na + 1)
    ac_row_ptr[1] = 1
    c_flag = zeros(Int, nc)
    for i = 1:na
        my_length = 0
        for k1 = ab_row_ptr[i]:ab_row_ptr[i + 1] - 1
            @inbounds b = ab_col_idx[k1]
            for k2 = bc_row_ptr[b]:bc_row_ptr[b + 1] - 1
                @inbounds c = bc_col_idx[k2]
                if (c_flag[c] < i)
                    @inbounds c_flag[c] = i
                    my_length += 1
                end
            end
        end
        ac_row_ptr[i + 1] = ac_row_ptr[i] + my_length
    end

    # Phase 2: Allocate column array
    nnz = ac_row_ptr[na + 1] - 1
    ac_col_idx = zeros(Int, nnz)

    # Phase 3: Fill entries
    fill!(c_flag, 0)
    for i = 1:na
        my_length = 0
        for k1 = ab_row_ptr[i]:ab_row_ptr[i + 1] - 1
            @inbounds b = ab_col_idx[k1]
            for k2 = bc_row_ptr[b]:bc_row_ptr[b + 1] - 1
                @inbounds c = bc_col_idx[k2]
                if (c_flag[c] < i)
                    @inbounds c_flag[c] = i
                    @inbounds ac_col_idx[ac_row_ptr[i] + my_length] = c
                    my_length += 1
                end
            end
        end
    end

    # Phase 4: Sort each row
    Threads.@threads for i = 1:na
        sort!(view(ac_col_idx, ac_row_ptr[i]:ac_row_ptr[i+1]-1))
    end

    return ac_row_ptr, ac_col_idx
end

"""
Greedy graph coloring algorithm (distance-1)
Returns: (color_to_nodes, node_to_color)
"""
function color_graph(n::Int, row_ptr::Vector{Int}, col_idx::Vector{Int})
    node_to_color = zeros(Int, n)
    max_color = 0

    for i = 1:n
        # Find colors used by neighbors
        neighbor_colors = Set{Int}()
        for k = row_ptr[i]:row_ptr[i + 1] - 1
            j = col_idx[k]
            if node_to_color[j] > 0
                push!(neighbor_colors, node_to_color[j])
            end
        end

        # Assign smallest available color
        color = 1
        while color in neighbor_colors
            color += 1
        end
        node_to_color[i] = color
        max_color = max(max_color, color)
    end

    # Build color_to_nodes
    color_to_nodes = [Int[] for _ = 1:max_color]
    for i = 1:n
        push!(color_to_nodes[node_to_color[i]], i)
    end

    return color_to_nodes, node_to_color
end

# ============================================================================
# Mesh Connectivity
# ============================================================================

"""
Build mesh connectivity graphs and coloring
Returns: (n2n_row_ptr, n2n_col_idx, e2e_colors)
"""
function build_mesh_connectivity(mesh::Mesh)
    nel = size(mesh.cell_to_node, 1)
    nnodes = length(mesh.vertex_x)

    # Build element-to-node graph (e2n)
    e2n_row_ptr = collect(1:4:4*nel+1)
    e2n_col_idx = collect(vec(mesh.cell_to_node'))

    # Transpose to get node-to-element (n2e)
    n2e_row_ptr, n2e_col_idx = transpose_graph(nel, nnodes, e2n_row_ptr, e2n_col_idx)

    # Combine to get node-to-node (n2n)
    n2n_row_ptr, n2n_col_idx = combine_graphs(nnodes, nel, nnodes,
                                               n2e_row_ptr, n2e_col_idx,
                                               e2n_row_ptr, e2n_col_idx)

    # Build element-to-element (e2e) for coloring
    e2e_row_ptr, e2e_col_idx = combine_graphs(nel, nnodes, nel,
                                               e2n_row_ptr, e2n_col_idx,
                                               n2e_row_ptr, n2e_col_idx)

    # Color the element graph
    e2e_colors, _ = color_graph(nel, e2e_row_ptr, e2e_col_idx)

    return n2n_row_ptr, n2n_col_idx, e2e_colors
end

# ============================================================================
# FEM Assembly
# ============================================================================

"""
Assemble element stiffness matrix and RHS for scaled Laplacian
Templated on element type parameters for compile-time optimization
Uses workspace to avoid allocations
"""
@inline function assemble_element!(elem::Q1Element{Dim, NNodes},
                          workspace::ElementWorkspace{Float64, Dim, NNodes},
                          x::Vector{Float64}, y::Vector{Float64},
                          ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    # Use workspace arrays (no allocation!)
    Ke = workspace.Ke
    fe = workspace.fe
    N = workspace.N
    dN_ref = workspace.dN_ref
    dN_phys = workspace.dN_phys
    J = workspace.J
    invJ = workspace.invJ

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    # Loop over quadrature points
    @fastmath @inbounds for (qp, (ξ, η)) in enumerate(elem.qpts)
        w = elem.qwts[qp]

        # Shape functions and gradients in reference coords (in-place)
        shape_functions!(N, ξ, η)
        shape_gradients!(dN_ref, ξ, η)

        # Compute Jacobian (manually unrolled for better performance)
        J[1, 1] = dN_ref[1, 1] * x[1] + dN_ref[2, 1] * x[2] + dN_ref[3, 1] * x[3] + dN_ref[4, 1] * x[4]
        J[1, 2] = dN_ref[1, 2] * x[1] + dN_ref[2, 2] * x[2] + dN_ref[3, 2] * x[3] + dN_ref[4, 2] * x[4]
        J[2, 1] = dN_ref[1, 1] * y[1] + dN_ref[2, 1] * y[2] + dN_ref[3, 1] * y[3] + dN_ref[4, 1] * y[4]
        J[2, 2] = dN_ref[1, 2] * y[1] + dN_ref[2, 2] * y[2] + dN_ref[3, 2] * y[3] + dN_ref[4, 2] * y[4]

        detJ = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]

        # Compute inverse Jacobian in-place
        invdetJ = 1.0 / detJ
        invJ[1, 1] =  J[2, 2] * invdetJ
        invJ[1, 2] = -J[1, 2] * invdetJ
        invJ[2, 1] = -J[2, 1] * invdetJ
        invJ[2, 2] =  J[1, 1] * invdetJ

        # Transform gradients to physical coordinates: dN/dx = dN/dξ * invJ
        @simd for i = 1:NNodes
            dN_phys[i, 1] = dN_ref[i, 1] * invJ[1, 1] + dN_ref[i, 2] * invJ[2, 1]
            dN_phys[i, 2] = dN_ref[i, 1] * invJ[1, 2] + dN_ref[i, 2] * invJ[2, 2]
        end

        # Physical coordinates at quadrature point
        xq = N[1] * x[1] + N[2] * x[2] + N[3] * x[3] + N[4] * x[4]
        yq = N[1] * y[1] + N[2] * y[2] + N[3] * y[3] + N[4] * y[4]

        # Evaluate coefficients
        ax_val = ax_func(xq, yq, 0.0)
        ay_val = ay_func(xq, yq, 0.0)
        f_val = f_func(xq, yq, 0.0)

        # Pre-compute common factors
        detJ_w = detJ * w
        ax_detJ_w = ax_val * detJ_w
        ay_detJ_w = ay_val * detJ_w
        f_detJ_w = f_val * detJ_w

        # Assemble stiffness matrix: Ke += (ax * dN/dx ⊗ dN/dx + ay * dN/dy ⊗ dN/dy) * detJ * w
        for i = 1:NNodes
            fe[i] += N[i] * f_detJ_w
            alphax_dNx_i = ax_detJ_w * dN_phys[i, 1]
            alphay_dNy_i = ay_detJ_w * dN_phys[i, 2]
            @simd for j = 1:NNodes
                Ke[i, j] += (alphax_dNx_i * dN_phys[j, 1] +
                            alphay_dNy_i * dN_phys[j, 2])
            end
        end
    end
end

# ============================================================================
# Sparse Matrix Utilities
# ============================================================================

"""
Find position of column j in row i of CSR matrix using binary search
Assumes col_idx is sorted within each row
"""
function find_matrix_position(row_ptr::Vector{Int}, col_idx::Vector{Int}, i::Int, j::Int)
    @inbounds left = row_ptr[i]
    @inbounds right = row_ptr[i + 1] - 1

    while left <= right
        mid = (left + right) >> 1
        @inbounds col_mid = col_idx[mid]
        if col_mid == j
            return mid
        elseif col_mid < j
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return -1
end

# ============================================================================
# Boundary Conditions
# ============================================================================

"""
Apply homogeneous Dirichlet BCs and build reduced system
"""
function apply_boundary_conditions(n2n_row_ptr::Vector{Int}, n2n_col_idx::Vector{Int},
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
        @inbounds gi = free_to_global[i]
        @inbounds reduced_rhs[i] = rhs[gi]
        for k = n2n_row_ptr[gi]:n2n_row_ptr[gi + 1] - 1
            gj = n2n_col_idx[k]
            if global_to_free[gj] != -1
                push!(reduced_col_idx, global_to_free[gj])
                push!(reduced_values, mat_values[k])
            end
        end
        @inbounds reduced_row_ptr[i + 1] = length(reduced_col_idx) + 1
    end

    return reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global
end

# ============================================================================
# Output
# ============================================================================

"""
Write solution to simple text format (coordinates + solution value)
"""
function write_solution(filename::String, mesh::Mesh, solution::Vector{Float64})
    open(filename, "w") do io
        println(io, "# x y u")
        for i = 1:length(solution)
            @printf(io, "%.6f %.6f %.6e\n", mesh.vertex_x[i], mesh.vertex_y[i], solution[i])
        end
    end
end

# ============================================================================
# MFEM Local Sparsity Pattern Construction
# ============================================================================

"""
    build_mfem_local_sparsity!(workspace, ratio_x, ratio_y) -> (nfree, nboundary, nnz_ii, nnz_b)

Build DOF mappings and sparsity patterns for MFEM static condensation of a single element.

This function constructs:
1. DOF mappings separating interior (free) and boundary nodes
2. K_ii sparsity pattern: interior-interior coupling (9-point stencil)
3. K_b sparsity pattern: boundary-all coupling

The fine grid has (ratio_x+1) × (ratio_y+1) nodes.

# Arguments
- `workspace`: Pre-allocated workspace containing arrays to fill
- `ratio_x::Int`: Refinement ratio in x-direction
- `ratio_y::Int`: Refinement ratio in y-direction

# Returns
- `nfree::Int`: Number of interior (free) DOFs
- `nboundary::Int`: Number of boundary DOFs
- `nnz_ii::Int`: Number of nonzeros in K_ii
- `nnz_b::Int`: Number of nonzeros in K_b

# Notes
- Fills workspace arrays in-place (no allocations)
- For square fine grids, use same value for ratio_x and ratio_y
- Node numbering: nodeID = ix + iy * (ratio_x + 1) + 1 (1-based)
"""
function build_mfem_local_sparsity!(workspace, ratio_x::Int, ratio_y::Int)
    # Extract DOF mapping arrays from workspace
    globalToFree = workspace.globalToFree
    freeToGlobal = workspace.freeToGlobal
    globalToBoundary = workspace.globalToBoundary
    boundaryToGlobal = workspace.boundaryToGlobal

    # Reset DOF mappings
    fill!(globalToFree, -1)
    fill!(globalToBoundary, -1)

    # Build DOF mappings: separate interior and boundary nodes
    nfree = 0
    nboundary = 0
    for iy = 0:ratio_y
        for ix = 0:ratio_x
            nodeID = ix + iy * (ratio_x + 1) + 1
            # Interior nodes
            if ix > 0 && ix < ratio_x && iy > 0 && iy < ratio_y
                nfree = nfree + 1
                @inbounds freeToGlobal[nfree] = nodeID
                @inbounds globalToFree[nodeID] = nfree
            else
                # Boundary nodes
                nboundary = nboundary + 1
                @inbounds boundaryToGlobal[nboundary] = nodeID
                @inbounds globalToBoundary[nodeID] = nboundary
            end
        end
    end

    # Extract K_ii arrays from workspace
    colptr_ii = workspace.colptr_ii
    rowidx_ii = workspace.rowidx_ii
    nzval_ii = workspace.nzval_ii

    # Reset K_ii values (sparsity pattern will be rebuilt)
    fill!(nzval_ii, 0.0)

    # Build K_ii sparsity pattern (interior-interior coupling)
    # Note: Building row-wise (CSR-style) but stored as CSC due to symmetry
    colptr_ii[1] = 1

    for i = 1:nfree
        @inbounds iGlobal = freeToGlobal[i]
        ix = (iGlobal - 1) % (ratio_x + 1)
        iy = (iGlobal - 1) ÷ (ratio_x + 1)

        # Count interior neighbors in 9-point stencil
        count = 1  # Diagonal
        hasWest = (ix > 1)
        hasEast = (ix < ratio_x - 1)
        hasSouth = (iy > 1)
        hasNorth = (iy < ratio_y - 1)

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
        ix = (iGlobal - 1) % (ratio_x + 1)
        iy = (iGlobal - 1) ÷ (ratio_x + 1)
        @inbounds offset = colptr_ii[i]

        # Add interior neighbors
        # South row
        if iy > 1
            if ix > 1
                jGlobal = iGlobal - 1 - (ratio_x + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
            jGlobal = iGlobal - (ratio_x + 1)
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
            if ix < ratio_x - 1
                jGlobal = iGlobal + 1 - (ratio_x + 1)
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
        if ix < ratio_x - 1
            jGlobal = iGlobal + 1
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
        end
        # North row
        if iy < ratio_y - 1
            if ix > 1
                jGlobal = iGlobal - 1 + (ratio_x + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
            jGlobal = iGlobal + (ratio_x + 1)
            @inbounds jFree = globalToFree[jGlobal]
            if jFree != -1
                @inbounds rowidx_ii[offset] = jFree
                offset += 1
            end
            if ix < ratio_x - 1
                jGlobal = iGlobal + 1 + (ratio_x + 1)
                @inbounds jFree = globalToFree[jGlobal]
                if jFree != -1
                    @inbounds rowidx_ii[offset] = jFree
                    offset += 1
                end
            end
        end
    end

    # Extract K_b arrays from workspace
    matRowPtr_b = workspace.matRowPtr_b
    matColIdx_b = workspace.matColIdx_b
    matValues_b = workspace.matValues_b

    # Reset K_b values (sparsity pattern will be rebuilt)
    fill!(matValues_b, 0.0)

    # Build K_b sparsity pattern (boundary-all coupling)
    @inbounds matRowPtr_b[1] = 1

    for i = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[i]
        ix = (iGlobal - 1) % (ratio_x + 1)
        iy = (iGlobal - 1) ÷ (ratio_x + 1)

        # Count all neighbors (boundary and interior) in global numbering
        count = 1  # Diagonal
        hasWest = (ix > 0)
        hasEast = (ix < ratio_x)

        count += hasWest + hasEast
        if iy > 0
            count += 1 + hasWest + hasEast
        end
        if iy < ratio_y
            count += 1 + hasWest + hasEast
        end

        @inbounds matRowPtr_b[i + 1] = matRowPtr_b[i] + count
    end

    @inbounds nnz_b = matRowPtr_b[nboundary + 1] - 1

    # Fill K_b column indices (in global node numbering)
    for i = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[i]
        ix = (iGlobal - 1) % (ratio_x + 1)
        iy = (iGlobal - 1) ÷ (ratio_x + 1)
        @inbounds offset = matRowPtr_b[i]

        # Add all neighbors
        # South row
        if iy > 0
            if ix > 0
                @inbounds matColIdx_b[offset] = iGlobal - 1 - (ratio_x + 1)
                offset += 1
            end
            @inbounds matColIdx_b[offset] = iGlobal - (ratio_x + 1)
            offset += 1
            if ix < ratio_x
                @inbounds matColIdx_b[offset] = iGlobal + 1 - (ratio_x + 1)
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
        if ix < ratio_x
            @inbounds matColIdx_b[offset] = iGlobal + 1
            offset += 1
        end
        # North row
        if iy < ratio_y
            if ix > 0
                @inbounds matColIdx_b[offset] = iGlobal - 1 + (ratio_x + 1)
                offset += 1
            end
            @inbounds matColIdx_b[offset] = iGlobal + (ratio_x + 1)
            offset += 1
            if ix < ratio_x
                @inbounds matColIdx_b[offset] = iGlobal + 1 + (ratio_x + 1)
                offset += 1
            end
        end
    end

    return nfree, nboundary, nnz_ii, nnz_b
end

end # module FEMBase
