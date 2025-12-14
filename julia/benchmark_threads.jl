#!/usr/bin/env julia
#
# Benchmark openmp_assembly.jl with varying thread counts
#

include("openmp_assembly.jl")

using Printf

function run_benchmark(nthreads::Int)
    println()
    println("="^70)
    println("BENCHMARK: $nthreads threads")
    println("="^70)
    println()

    # Problem setup
    ax(x, y, z) = 1.0 + 0.5 * sin(2π * x) * cos(2π * y)
    ay(x, y, z) = 1.0 + 0.5 * cos(2π * x) * sin(2π * y)
    f(x, y, z) = 2.0 * π^2 * sin(π * x) * sin(π * y)

    # Mesh generation
    nx, ny = 512, 512
    t0 = time()
    mesh = generate_mesh(nx, ny)
    mesh_time = time() - t0

    # Connectivity
    t0 = time()
    n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
    conn_time = time() - t0

    # Assembly
    elem = Q1Element()
    t0 = time()
    mat_values, rhs, color_times, lookup_time = assemble_system(
        mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f
    )
    assembly_time = time() - t0

    # Boundary conditions
    t0 = time()
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global =
        apply_boundary_conditions(n2n_row_ptr, n2n_col_idx, mat_values, rhs, mesh.boundary_nodes)
    bc_time = time() - t0

    # Solve
    nfree = length(free_to_global)
    K = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)
    t0 = time()
    sol_free = K \ reduced_rhs
    solve_time = time() - t0

    # Report
    total_time = mesh_time + conn_time + assembly_time + bc_time + solve_time

    println("Performance Summary:")
    println("  Mesh generation:     ", @sprintf("%8.2f ms", mesh_time * 1000))
    println("  Connectivity/color:  ", @sprintf("%8.2f ms", conn_time * 1000))
    println("  Assembly (total):    ", @sprintf("%8.2f ms", assembly_time * 1000))
    println("    Lookup table:      ", @sprintf("%8.2f ms", lookup_time * 1000))
    println("    Element assembly:  ", @sprintf("%8.2f ms", (assembly_time - lookup_time) * 1000))
    println("  Boundary conditions: ", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Solve:               ", @sprintf("%8.2f ms", solve_time * 1000))
    println("  " * "-"^68)
    println("  Total:               ", @sprintf("%8.2f ms", total_time * 1000))
    println()

    println("  Color breakdown:")
    for ic = 1:length(e2e_colors)
        println("    Color $ic (", @sprintf("%6d", length(e2e_colors[ic])), " elems): ",
                @sprintf("%8.2f ms", color_times[ic] * 1000))
    end
    println()

    return Dict(
        "nthreads" => nthreads,
        "mesh_time" => mesh_time * 1000,
        "conn_time" => conn_time * 1000,
        "assembly_time" => assembly_time * 1000,
        "lookup_time" => lookup_time * 1000,
        "element_time" => (assembly_time - lookup_time) * 1000,
        "bc_time" => bc_time * 1000,
        "solve_time" => solve_time * 1000,
        "total_time" => total_time * 1000,
        "color_times" => color_times * 1000
    )
end

function main()
    println("="^70)
    println("Multi-threaded Performance Benchmark")
    println("="^70)
    println()
    println("Problem: 512x512 Q1 mesh, scaled Laplacian")
    println()

    # Get available threads
    max_threads = Threads.nthreads()
    println("Available threads: $max_threads")

    if max_threads < 8
        println("WARNING: Running with only $max_threads threads.")
        println("         For best results, run with: JULIA_NUM_THREADS=8 julia benchmark_threads.jl")
    end
    println()

    # Run benchmarks
    thread_counts = [1, 2, 4, 8]
    thread_counts = filter(n -> n <= max_threads, thread_counts)

    results = []

    # Warmup run (compile everything)
    println("Warmup run...")
    run_benchmark(max_threads)

    # Actual benchmarks
    for nt in thread_counts
        # Note: Can't change thread count at runtime in Julia
        # This will use all available threads
        result = run_benchmark(max_threads)
        push!(results, result)
    end

    # Summary table
    println()
    println("="^70)
    println("SUMMARY TABLE")
    println("="^70)
    println()
    println("Note: Julia cannot change thread count at runtime.")
    println("      All runs use $max_threads threads.")
    println("      To test different thread counts, restart Julia with:")
    println("      JULIA_NUM_THREADS=N julia benchmark_threads.jl")
    println()
    println(@sprintf("%-10s %10s %10s %10s %10s %10s",
                     "Threads", "Mesh", "Connect", "Assembly", "Solve", "Total"))
    println("-"^70)

    for r in results
        println(@sprintf("%-10d %9.1f %10.1f %10.1f %10.1f %10.1f",
                        r["nthreads"],
                        r["mesh_time"],
                        r["conn_time"],
                        r["assembly_time"],
                        r["solve_time"],
                        r["total_time"]))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
