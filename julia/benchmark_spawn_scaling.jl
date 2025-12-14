#!/usr/bin/env julia
#
# Test @spawn vs @threads across different thread counts
#

include("openmp_assembly.jl")
using Printf
using Statistics

function quick_benchmark(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, method, nchunks)
    # 3 trials
    times = []
    for _ = 1:3
        if method == :threads
            t0 = time()
            assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
            push!(times, (time() - t0) * 1000)
        else
            t0 = time()
            assemble_system_spawn(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f;
                                 nchunks_per_color=nchunks)
            push!(times, (time() - t0) * 1000)
        end
    end
    return median(times)
end

# Setup
ax(x, y, z) = 1.0 + 0.5 * sin(2π * x) * cos(2π * y)
ay(x, y, z) = 1.0 + 0.5 * cos(2π * x) * sin(2π * y)
f(x, y, z) = 2.0 * π^2 * sin(π * x) * sin(π * y)

nx, ny = 512, 512
mesh = generate_mesh(nx, ny)
n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
elem = Q1Element()

nthreads = Threads.nthreads()

println("Thread count: $nthreads")
println("Testing assembly methods...")
println()

# Warmup
assemble_system(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
assemble_system_spawn(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)

# Test with optimal chunk count (1× threads)
optimal_chunks = nthreads

threads_time = quick_benchmark(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, :threads, 0)
spawn_time = quick_benchmark(mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f, :spawn, optimal_chunks)

speedup = threads_time / spawn_time

println("Results (median of 3 trials):")
println(@sprintf("  @threads:              %8.2f ms", threads_time))
println(@sprintf("  @spawn (%d chunks):    %8.2f ms", optimal_chunks, spawn_time))
println(@sprintf("  Speedup:               %8.3fx (%.1f%% %s)",
                 speedup, abs(speedup - 1.0) * 100,
                 speedup > 1.0 ? "faster" : "slower"))
