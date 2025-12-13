#!/usr/bin/env julia
#
# @file openmp_assembly.jl
# @brief Julia version of OpenMP FEM assembly example
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

# ============================================================================
# Q1 Finite Element (2D Bilinear Quadrilateral)
# ============================================================================

"""
Q1 shape functions and gradients on reference element [-1,1]²
"""
struct Q1Element
    # Gauss quadrature points and weights for 2x2 rule
    qpts::Vector{Tuple{Float64, Float64}}
    qwts::Vector{Float64}

    function Q1Element()
        # 2x2 Gauss quadrature
        gp = 1.0 / sqrt(3.0)
        qpts = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
        qwts = [1.0, 1.0, 1.0, 1.0]
        new(qpts, qwts)
    end
end

"""
Evaluate Q1 shape functions at (ξ, η)
Returns: Vector of 4 shape function values [N1, N2, N3, N4]
Node ordering: (ξ,η) = [(-1,-1), (1,-1), (1,1), (-1,1)]
"""
function shape_functions(ξ::Float64, η::Float64)
    N = zeros(4)
    N[1] = 0.25 * (1 - ξ) * (1 - η)
    N[2] = 0.25 * (1 + ξ) * (1 - η)
    N[3] = 0.25 * (1 + ξ) * (1 + η)
    N[4] = 0.25 * (1 - ξ) * (1 + η)
    return N
end

"""
Evaluate Q1 shape function gradients in reference coordinates
Returns: 4x2 matrix where row i contains [∂Ni/∂ξ, ∂Ni/∂η]
"""
function shape_gradients(ξ::Float64, η::Float64)
    dN = zeros(4, 2)
    # ∂N/∂ξ
    dN[1, 1] = -0.25 * (1 - η)
    dN[2, 1] =  0.25 * (1 - η)
    dN[3, 1] =  0.25 * (1 + η)
    dN[4, 1] = -0.25 * (1 + η)
    # ∂N/∂η
    dN[1, 2] = -0.25 * (1 - ξ)
    dN[2, 2] = -0.25 * (1 + ξ)
    dN[3, 2] =  0.25 * (1 + ξ)
    dN[4, 2] =  0.25 * (1 - ξ)
    return dN
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
Transpose a graph stored as CSR (row_ptr, col_idx)
"""
function transpose_graph(n::Int, m::Int, row_ptr::Vector{Int}, col_idx::Vector{Int})
    nnz = length(col_idx)

    # Count entries per column
    col_count = zeros(Int, m)
    for j in col_idx
        col_count[j] += 1
    end

    # Build new row_ptr
    new_row_ptr = zeros(Int, m + 1)
    new_row_ptr[1] = 1
    for i = 1:m
        new_row_ptr[i + 1] = new_row_ptr[i] + col_count[i]
    end

    # Fill new_col_idx
    new_col_idx = zeros(Int, nnz)
    col_offset = zeros(Int, m)

    for i = 1:n
        for k = row_ptr[i]:row_ptr[i + 1] - 1
            j = col_idx[k]
            pos = new_row_ptr[j] + col_offset[j]
            new_col_idx[pos] = i
            col_offset[j] += 1
        end
    end

    return new_row_ptr, new_col_idx
end

"""
Combine two graphs: A → B and B → C to get A → C
"""
function combine_graphs(na::Int, nb::Int, nc::Int,
                       ab_row_ptr::Vector{Int}, ab_col_idx::Vector{Int},
                       bc_row_ptr::Vector{Int}, bc_col_idx::Vector{Int})
    # Result: ac[i] = union of bc[ab[i][j]] for all j
    ac_sets = [Set{Int}() for _ = 1:na]

    for i = 1:na
        for k1 = ab_row_ptr[i]:ab_row_ptr[i + 1] - 1
            b = ab_col_idx[k1]
            for k2 = bc_row_ptr[b]:bc_row_ptr[b + 1] - 1
                c = bc_col_idx[k2]
                push!(ac_sets[i], c)
            end
        end
    end

    # Convert sets to CSR
    ac_row_ptr = zeros(Int, na + 1)
    ac_row_ptr[1] = 1
    for i = 1:na
        ac_row_ptr[i + 1] = ac_row_ptr[i] + length(ac_sets[i])
    end

    nnz = ac_row_ptr[na + 1] - 1
    ac_col_idx = zeros(Int, nnz)

    pos = 1
    for i = 1:na
        for c in sort(collect(ac_sets[i]))
            ac_col_idx[pos] = c
            pos += 1
        end
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
"""
function assemble_element!(Ke::Matrix{Float64}, fe::Vector{Float64},
                          elem::Q1Element,
                          x::Vector{Float64}, y::Vector{Float64},
                          ax_func::Function, ay_func::Function, f_func::Function)
    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    # Loop over quadrature points
    for (qp, (ξ, η)) in enumerate(elem.qpts)
        w = elem.qwts[qp]

        # Shape functions and gradients in reference coords
        N = shape_functions(ξ, η)
        dN_ref = shape_gradients(ξ, η)  # 4x2: [∂N/∂ξ, ∂N/∂η]

        # Compute Jacobian
        J = zeros(2, 2)
        for i = 1:4
            J[1, 1] += dN_ref[i, 1] * x[i]  # ∂x/∂ξ
            J[1, 2] += dN_ref[i, 2] * x[i]  # ∂x/∂η
            J[2, 1] += dN_ref[i, 1] * y[i]  # ∂y/∂ξ
            J[2, 2] += dN_ref[i, 2] * y[i]  # ∂y/∂η
        end

        detJ = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
        invJ = [J[2, 2] -J[1, 2]; -J[2, 1] J[1, 1]] / detJ

        # Transform gradients to physical coordinates: dN/dx = invJ^T * dN/dξ
        dN_phys = dN_ref * invJ'  # 4x2: [∂N/∂x, ∂N/∂y]

        # Physical coordinates at quadrature point
        xq = sum(N[i] * x[i] for i = 1:4)
        yq = sum(N[i] * y[i] for i = 1:4)

        # Evaluate coefficients
        ax_val = ax_func(xq, yq, 0.0)
        ay_val = ay_func(xq, yq, 0.0)
        f_val = f_func(xq, yq, 0.0)

        # Assemble stiffness matrix: Ke += (ax * dN/dx ⊗ dN/dx + ay * dN/dy ⊗ dN/dy) * detJ * w
        for i = 1:4
            for j = 1:4
                Ke[i, j] += (ax_val * dN_phys[i, 1] * dN_phys[j, 1] +
                            ay_val * dN_phys[i, 2] * dN_phys[j, 2]) * detJ * w
            end
        end

        # Assemble RHS: fe += N * f * detJ * w
        for i = 1:4
            fe[i] += N[i] * f_val * detJ * w
        end
    end
end

"""
Find position of column j in row i of CSR matrix
"""
function find_matrix_position(row_ptr::Vector{Int}, col_idx::Vector{Int}, i::Int, j::Int)
    for k = row_ptr[i]:row_ptr[i + 1] - 1
        if col_idx[k] == j
            return k
        end
    end
    return -1
end

"""
Assemble global system with thread-parallel coloring
"""
function assemble_system(mesh::Mesh, elem::Q1Element,
                        n2n_row_ptr::Vector{Int}, n2n_col_idx::Vector{Int},
                        e2e_colors::Vector{Vector{Int}},
                        ax_func::Function, ay_func::Function, f_func::Function)
    nnodes = length(mesh.vertex_x)
    nnz = length(n2n_col_idx)

    # Initialize global arrays
    mat_values = zeros(nnz)
    rhs = zeros(nnodes)

    # Thread-local element matrices
    num_threads = Threads.nthreads()
    Ke_local = [zeros(4, 4) for _ = 1:num_threads]
    fe_local = [zeros(4) for _ = 1:num_threads]

    # Loop over colors
    for (ic, elements) in enumerate(e2e_colors)
        # All elements in this color can be assembled in parallel (no conflicts)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            Ke = Ke_local[tid]
            fe = fe_local[tid]

            # Get element nodes
            nodes = mesh.cell_to_node[iel, :]
            x = mesh.vertex_x[nodes]
            y = mesh.vertex_y[nodes]

            # Assemble element
            assemble_element!(Ke, fe, elem, x, y, ax_func, ay_func, f_func)

            # Scatter to global arrays (no race condition within same color)
            for i = 1:4
                gi = nodes[i]
                rhs[gi] += fe[i]

                for j = 1:4
                    gj = nodes[j]
                    k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                    if k > 0
                        mat_values[k] += Ke[i, j]
                    end
                end
            end
        end
    end

    return mat_values, rhs
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
        gi = free_to_global[i]
        reduced_rhs[i] = rhs[gi]

        for k = n2n_row_ptr[gi]:n2n_row_ptr[gi + 1] - 1
            gj = n2n_col_idx[k]
            if global_to_free[gj] != -1
                push!(reduced_col_idx, global_to_free[gj])
                push!(reduced_values, mat_values[k])
            end
        end
        reduced_row_ptr[i + 1] = length(reduced_col_idx) + 1
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

    # Material coefficients and forcing term
    ax(x, y, z) = 1.0
    ay(x, y, z) = 1.0
    function f(x, y, z)
        # Manufactured solution: u = sin(π*x) * sin(π*y)
        # => -Δu = 2π² * sin(π*x) * sin(π*y)
        return 2.0 * π^2 * sin(π * x) * sin(π * y)
    end

    # ========================================================================
    # Mesh generation
    # ========================================================================

    nx, ny = 64, 64
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
    # Assembly
    # ========================================================================

    println("="^70)
    println("Starting Assembly")
    println("="^70)

    elem = Q1Element()

    t0 = time()
    mat_values, rhs = assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors,
                                     ax, ay, f)
    assembly_time = time() - t0

    println("Assembly complete")
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
