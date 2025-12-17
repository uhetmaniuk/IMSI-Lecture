#!/usr/bin/env julia
#
# @file openmp_assembly.jl
# @brief Julia version of OpenMP FEM assembly example (Optimized)
#
# Solves: -∇·(α∇u) = f on a 2D rectangular domain
# with homogeneous Dirichlet boundary conditions
#
# Features:
# - 2D Q1 (bilinear quadrilateral) finite elements
# - Varying coefficients ax, ay, f
# - Thread-parallel assembly with graph coloring
# - Built-in sparse solver (\)
#

using SparseArrays
using LinearAlgebra
using Printf
using PrecompileTools: @compile_workload

# ============================================================================
# Q1 Finite Element (2D Bilinear Quadrilateral)
# ============================================================================

"""
Q1 shape functions and gradients on reference element [-1,1]²
Type parameters Dim and NNodes are compile-time constants
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

# Accessor methods to maintain compatibility
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
    # Use maxthreadid() to get the maximum possible thread ID (accounts for multiple thread pools)
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

    # Reduce counts from all threads (only up to max_tid that was actually used)
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
Workspace for element assembly (reused across elements to avoid allocations)
"""
struct ElementWorkspace{T<:AbstractFloat, Dim, NNodes}
    N::Vector{T}
    dN_ref::Matrix{T}
    dN_phys::Matrix{T}
    J::Matrix{T}
    invJ::Matrix{T}

    function ElementWorkspace{T, Dim, NNodes}() where {T<:AbstractFloat, Dim, NNodes}
        new{T, Dim, NNodes}(
            zeros(T, NNodes),
            zeros(T, NNodes, Dim),
            zeros(T, NNodes, Dim),
            zeros(T, Dim, Dim),
            zeros(T, Dim, Dim)
        )
    end
end

"""
Assemble element stiffness matrix and RHS for scaled Laplacian
Templated on element type parameters for compile-time optimization
Uses workspace to avoid allocations
"""
@inline function assemble_element!(Ke::Matrix{Float64}, fe::Vector{Float64},
                          elem::Q1Element{Dim, NNodes},
                          workspace::ElementWorkspace{Float64, Dim, NNodes},
                          x::Vector{Float64}, y::Vector{Float64},
                          ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    # Use workspace arrays (no allocation!)
    N = workspace.N
    dN_ref = workspace.dN_ref
    dN_phys = workspace.dN_phys
    J = workspace.J
    invJ = workspace.invJ

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
            dNx_i = dN_phys[i, 1]
            dNy_i = dN_phys[i, 2]
            fe[i] += N[i] * f_detJ_w

            @simd for j = 1:NNodes
                Ke[i, j] += (ax_detJ_w * dNx_i * dN_phys[j, 1] +
                            ay_detJ_w * dNy_i * dN_phys[j, 2])
            end
        end
    end
end

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

"""
Assemble global system with thread-parallel coloring
Returns: mat_values, rhs, color_times, thread_stats
Templated on element type parameters for compile-time optimization
"""
function assemble_system(mesh::Mesh, elem::Q1Element{Dim, NNodes},
                        n2n_row_ptr::Vector{Int}, n2n_col_idx::Vector{Int},
                        e2e_colors::Vector{Vector{Int}},
                        ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    nnodes = length(mesh.vertex_x)
    nnz = length(n2n_col_idx)

    # Initialize global arrays
    mat_values = zeros(nnz)
    rhs = zeros(nnodes)

    # Thread-local element matrices (using type parameters for compile-time sizes)
    max_tid = Threads.maxthreadid()
    Ke_local = [zeros(NNodes, NNodes) for _ = 1:max_tid]
    fe_local = [zeros(NNodes) for _ = 1:max_tid]
    workspace_local = [ElementWorkspace{Float64, Dim, NNodes}() for _ = 1:max_tid]

    # Track time per color
    color_times = zeros(length(e2e_colors))

    # Loop over colors
    for (ic, elements) in enumerate(e2e_colors)
        t_color = time()

        # All elements in this color can be assembled in parallel (no conflicts)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            @inbounds Ke = Ke_local[tid]
            @inbounds fe = fe_local[tid]

            # Get element nodes
            @inbounds nodes = mesh.cell_to_node[iel, :]
            @inbounds x = mesh.vertex_x[nodes]
            @inbounds y = mesh.vertex_y[nodes]

            # Assemble element
            @inbounds workspace = workspace_local[tid]
            assemble_element!(Ke, fe, elem, workspace, x, y, ax_func, ay_func, f_func)

            # Scatter to global arrays
            # Scatter to global arrays (no race condition within same color)
            for i = 1:NNodes
                @inbounds gi = nodes[i]
                @inbounds rhs[gi] += fe[i]
            end

            for i = 1:NNodes
                @inbounds gi = nodes[i]
                for j = 1:NNodes
                    @inbounds gj = nodes[j]
                    k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                    @inbounds mat_values[k] += Ke[i, j]
                end
            end

        end

        @inbounds color_times[ic] = time() - t_color
    end

    return mat_values, rhs, color_times
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
# Precompilation Workload
# ============================================================================

@compile_workload begin
    # Small warm-up problem to precompile hot functions
    # This eliminates JIT overhead from first assembly iteration

    # Generate tiny mesh (4x4 elements)
    mesh_warm = generate_mesh(4, 4)
    elem_warm = Q1Element()

    # Build connectivity
    n2n_row_ptr_warm, n2n_col_idx_warm, e2e_colors_warm = build_mesh_connectivity(mesh_warm)

    # Simple coefficient functions
    ax_warm(x, y, z) = 1.0
    ay_warm(x, y, z) = 1.0
    f_warm(x, y, z) = 1.0

    # Run assembly (this precompiles all kernels)
    mat_values_warm, rhs_warm, _ = assemble_system(
        mesh_warm, elem_warm,
        n2n_row_ptr_warm, n2n_col_idx_warm, e2e_colors_warm,
        ax_warm, ay_warm, f_warm
    )

    # Precompile boundary conditions
    reduced_row_ptr_warm, reduced_col_idx_warm, reduced_values_warm, reduced_rhs_warm, free_to_global_warm =
        apply_boundary_conditions(n2n_row_ptr_warm, n2n_col_idx_warm, mat_values_warm, rhs_warm, mesh_warm.boundary_nodes)

    # Precompile sparse solve
    K_warm = SparseMatrixCSC(length(free_to_global_warm), length(free_to_global_warm),
                              reduced_row_ptr_warm, reduced_col_idx_warm, reduced_values_warm)
    sol_warm = K_warm \ reduced_rhs_warm
end

# ============================================================================
# Main Program
# ============================================================================

function main()
    println("="^70)
    println("Julia FEM Assembly Example (Thread-Parallel with Coloring)")
    println("="^70)
    println()

    # Thread info
    println("Number of threads: ", Threads.nthreads())
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

    nx, ny = 256, 256
    println("Generating mesh: $nx x $ny Q1 elements")

    t0 = time()
    mesh = generate_mesh(nx, ny)
    mesh_time = time() - t0

    println("  Number of elements: ", size(mesh.cell_to_node, 1))
    println("  Number of nodes:    ", length(mesh.vertex_x))
    println("  Mesh generation:    ", @sprintf("%.2f ms", mesh_time * 1000))
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
    println("  Connectivity time: ", @sprintf("%.2f ms", conn_time * 1000))
    println()

    # ========================================================================
    # Assembly (with JIT warm-up)
    # ========================================================================

    println("="^70)
    println("Starting Assembly")
    println("="^70)

    elem = Q1Element()

    # Warm-up run to eliminate JIT compilation overhead
    println("Running JIT warm-up...")
    t_warmup = time()
    _ = assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
    warmup_time = time() - t_warmup
    println("  Warm-up time: ", @sprintf("%.2f ms", warmup_time * 1000))
    println()

    # Actual timed run (using binary search for matrix position finding)
    println("Using binary search for matrix position finding")
    t0 = time()
    mat_values, rhs, color_times = assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors,
                                                                  ax, ay, f)
    assembly_time = time() - t0

    println("Assembly complete")
    println("  Total assembly time: ", @sprintf("%.2f ms", assembly_time * 1000))
    println()
    println("  Time per color:")
    for ic = 1:length(e2e_colors)
        println("    Color $ic (", length(e2e_colors[ic]), " elements): ",
                @sprintf("%.2f ms", color_times[ic] * 1000))
    end
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
    # Output
    # ========================================================================

    println("Writing solution to file...")

    # Expand to full solution (BCs = 0)
    solution = zeros(length(mesh.vertex_x))
    for i = 1:nfree
        solution[free_to_global[i]] = sol_free[i]
    end

    write_solution("julia/solution.txt", mesh, solution)
    println("  Solution written to: julia/solution.txt")
    println()

    # ========================================================================
    # Performance summary
    # ========================================================================

    println("="^70)
    println("Performance Summary")
    println("="^70)
    println("  Mesh generation:    ", @sprintf("%8.2f ms", mesh_time * 1000))
    println("  Connectivity/color: ", @sprintf("%8.2f ms", conn_time * 1000))
    println("  Assembly:           ", @sprintf("%8.2f ms", assembly_time * 1000))
    println("  Boundary conditions:", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Solve:              ", @sprintf("%8.2f ms", solve_time * 1000))
    println("  " * "-"^68)
    total_time = mesh_time + conn_time + assembly_time + bc_time + solve_time
    println("  Total:              ", @sprintf("%8.2f ms", total_time * 1000))
    println()

    # Solution statistics
    minval = minimum(sol_free)
    maxval = maximum(sol_free)
    avgval = sum(sol_free) / length(sol_free)

    println("Solution statistics:")
    println("  Min value: ", @sprintf("%.6e", minval))
    println("  Max value: ", @sprintf("%.6e", maxval))
    println("  Avg value: ", @sprintf("%.6e", avgval))
    println()
end

# Run main program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
