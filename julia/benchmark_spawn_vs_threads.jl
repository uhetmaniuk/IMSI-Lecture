#!/usr/bin/env julia
#
# Benchmark @spawn vs @threads for FEM assembly
#

include("openmp_assembly.jl")

using Printf
using Statistics

function benchmark_assembly_method(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors,
                                   ax, ay, f, method::Symbol, nchunks::Int=0; ntrials::Int=5)
    times = Float64[]
    lookup_times = Float64[]

    for trial = 1:ntrials
        if method == :threads
            t0 = time()
            mat_values, rhs, color_times, lookup_time = assemble_system(
                mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f
            )
            elapsed = time() - t0
        elseif method == :spawn
            t0 = time()
            mat_values, rhs, color_times, lookup_time = assemble_system_spawn(
                mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f;
                nchunks_per_color=nchunks
            )
            elapsed = time() - t0
        else
            error("Unknown method: $method")
        end

        push!(times, elapsed * 1000)
        push!(lookup_times, lookup_time * 1000)
    end

    return median(times), median(lookup_times), times
end

function main()
    println("="^70)
    println("Benchmark: @spawn vs @threads Assembly")
    println("="^70)
    println()

    # Problem setup
    ax(x, y, z) = 1.0 + 0.5 * sin(2π * x) * cos(2π * y)
    ay(x, y, z) = 1.0 + 0.5 * cos(2π * x) * sin(2π * y)
    f(x, y, z) = 2.0 * π^2 * sin(π * x) * sin(π * y)

    # Mesh
    nx, ny = 512, 512
    println("Mesh: $nx × $ny Q1 elements")
    println("Threads: ", Threads.nthreads())
    println()

    mesh = generate_mesh(nx, ny)
    n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
    elem = Q1Element()

    println("Number of colors: ", length(e2e_colors))
    println("Elements per color: ", length(e2e_colors[1]))
    println()

    # Warmup
    println("Warmup runs...")
    assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
    assemble_system_spawn(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
    println()

    # Benchmark configurations
    nthreads = Threads.nthreads()
    chunk_configs = [
        nthreads,      # 1x threads (coarse-grained)
        nthreads * 2,  # 2x threads (default)
        nthreads * 4,  # 4x threads (fine-grained)
        nthreads * 8,  # 8x threads (very fine)
    ]

    println("="^70)
    println("RESULTS")
    println("="^70)
    println()

    # Baseline: @threads
    println("Baseline: @threads")
    threads_time, threads_lookup, threads_trials = benchmark_assembly_method(
        mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, :threads
    )
    println(@sprintf("  Median time:   %.2f ms", threads_time))
    println(@sprintf("  Lookup table:  %.2f ms", threads_lookup))
    println(@sprintf("  Std dev:       %.2f ms", std(threads_trials)))
    println(@sprintf("  Min/Max:       %.2f / %.2f ms", minimum(threads_trials), maximum(threads_trials)))
    println()

    # Test @spawn with different chunk sizes
    println("@spawn with varying chunk counts:")
    println("-"^70)
    println(@sprintf("%-20s %12s %12s %12s", "Configuration", "Time (ms)", "Speedup", "vs Baseline"))
    println("-"^70)

    best_speedup = 1.0
    best_config = "baseline"

    for nchunks in chunk_configs
        spawn_time, spawn_lookup, spawn_trials = benchmark_assembly_method(
            mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, :spawn, nchunks
        )

        speedup = threads_time / spawn_time
        config_name = "$nchunks chunks/color"

        println(@sprintf("%-20s %12.2f %12.3fx %11.1f%%",
                        config_name, spawn_time, speedup, (speedup - 1.0) * 100))

        if speedup > best_speedup
            best_speedup = speedup
            best_config = config_name
        end
    end

    println("-"^70)
    println()

    # Summary
    println("="^70)
    println("SUMMARY")
    println("="^70)
    println()
    println("Baseline (@threads):        ", @sprintf("%.2f ms", threads_time))

    if best_speedup > 1.0
        println("Best configuration:         ", best_config)
        println("Best speedup:               ", @sprintf("%.3fx (%.1f%% faster)",
                                                          best_speedup, (best_speedup - 1.0) * 100))
    else
        println("Result: @threads is faster than all @spawn configurations")
        println("Reason: Likely due to:")
        println("  - Task spawning overhead exceeds benefits")
        println("  - Memory bandwidth is the bottleneck, not scheduling")
        println("  - Uniform workload doesn't benefit from dynamic scheduling")
    end
    println()

    # Detailed breakdown for best @spawn config
    println("="^70)
    println("DETAILED COMPARISON (Best @spawn config)")
    println("="^70)
    println()

    best_nchunks = chunk_configs[argmax([threads_time / benchmark_assembly_method(
        mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, :spawn, nc
    )[1] for nc in chunk_configs])]

    # Run detailed comparison
    println("Running detailed trials...")
    println()

    # @threads detailed
    threads_detailed = []
    for i = 1:10
        mat_values, rhs, color_times, lookup_time = assemble_system(
            mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f
        )
        push!(threads_detailed, sum(color_times) * 1000)
    end

    # @spawn detailed
    spawn_detailed = []
    for i = 1:10
        mat_values, rhs, color_times, lookup_time = assemble_system_spawn(
            mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f;
            nchunks_per_color=best_nchunks
        )
        push!(spawn_detailed, sum(color_times) * 1000)
    end

    println("@threads (element assembly only):")
    println(@sprintf("  Mean:   %.2f ms", mean(threads_detailed)))
    println(@sprintf("  Median: %.2f ms", median(threads_detailed)))
    println(@sprintf("  Stddev: %.2f ms", std(threads_detailed)))
    println()

    println("@spawn with $best_nchunks chunks (element assembly only):")
    println(@sprintf("  Mean:   %.2f ms", mean(spawn_detailed)))
    println(@sprintf("  Median: %.2f ms", median(spawn_detailed)))
    println(@sprintf("  Stddev: %.2f ms", std(spawn_detailed)))
    println()

    improvement = (1.0 - median(spawn_detailed) / median(threads_detailed)) * 100
    if improvement > 0
        println(@sprintf("Assembly speedup: %.1f%% faster with @spawn", improvement))
    else
        println(@sprintf("Assembly slowdown: %.1f%% slower with @spawn", -improvement))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
