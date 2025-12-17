#!/usr/bin/env julia
#
# @file q1_assembly.jl
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

# Load shared FEM base kernels
include("FEMBase.jl")
using .FEMBase

# ============================================================================
# System Assembly
# ============================================================================

"""
Assemble global system with thread-parallel coloring
Returns: mat_values, rhs, color_times
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

    write_solution("solution.txt", mesh, solution)
    println("  Solution written to: solution.txt")
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
