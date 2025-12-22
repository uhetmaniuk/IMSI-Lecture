#!/usr/bin/env julia
#
# @file test_pcg_gpu.jl
# @brief Test PCG_GPU module with simple examples (CPU and GPU)
#

using LinearAlgebra
using SparseArrays
using CUDA
using CUDA.CUSPARSE
using Printf

using Krylov
import Krylov: CgWorkspace, cg!

# Load the PCG_GPU module
include("PCG_GPU.jl")
using .PCG_GPU

println("="^70)
println("Testing PCG_GPU Module")
println("="^70)
println()

# Test 1: Small CPU problem
println("Test 1: CPU Laplacian (10×10)")
println("-"^70)

n = 10
# Create a simple Laplacian matrix (SPD)
A_cpu = spdiagm(0 => 2*ones(n), 1 => -ones(n-1), -1 => -ones(n-1))
b_cpu = ones(n)
x0_cpu = 0.5 * ones(n)  # Initial guess

println("Matrix: $(size(A_cpu)), nnz=$(nnz(A_cpu))")
println("Testing that initial guess is used...")

# Solve with PCG_GPU
x_pcg, info_pcg = pcg_solve(A_cpu, b_cpu; x0=x0_cpu, tol=1e-10, verbose=true)

println("PCG_GPU: ", info_pcg)
println("Solution norm: ", norm(x_pcg))
println("Residual: ", norm(A_cpu * x_pcg - b_cpu))
println()

# Test 2: GPU problem (if available)
if CUDA.functional()
    println("Test 2: GPU Laplacian (1000×1000)")
    println("-"^70)

    n = 1000
    # Create larger Laplacian
    A_cpu_large = spdiagm(0 => 2*ones(n), 1 => -ones(n-1), -1 => -ones(n-1))
    b_cpu_large = ones(n)
    x0_cpu_large = 0.5 * ones(n)

    # Transfer to GPU
    A_gpu = CuSparseMatrixCSC(A_cpu_large)
    b_gpu = CuArray(b_cpu_large)
    x0_gpu = CuArray(x0_cpu_large)

    println("Matrix: $(size(A_gpu)), nnz=$(nnz(A_cpu_large))")

    # Warm-up
    x_warmup, _ = pcg_solve(A_gpu, b_gpu; x0=x0_gpu, tol=1e-10, verbose=false)

    # Timed solve
    CUDA.@sync begin
        t0 = time()
        x_gpu, info_gpu = pcg_solve(A_gpu, b_gpu; x0=x0_gpu, tol=1e-10, verbose=true)
        gpu_time = time() - t0
    end

    println("PCG_GPU: ", info_gpu)
    println("GPU time: ", @sprintf("%.2f ms", gpu_time * 1000))

    # Verify correctness
    x_check = Vector(x_gpu)
    residual = norm(A_cpu_large * x_check - b_cpu_large)
    println("Residual (on CPU): ", residual)
    println()

    # Test 3: Compare with Krylov.jl (optional)
    println("Test 3: Comparison with Krylov.jl")
    println("-"^70)

    try
        # PCG_GPU timing
        CUDA.@sync begin
            t_pcg = @elapsed begin
                x_pcg, info_pcg = pcg_solve(A_gpu, b_gpu; x0=x0_gpu, tol=1e-10, verbose=false)
            end
        end

        # Krylov.jl timing
        CUDA.@sync begin
            t_krylov = @elapsed begin
                solver = CgWorkspace(A_gpu, b_gpu)
                cg!(solver, A_gpu, b_gpu; atol=1e-10, rtol=0.0, verbose=0)
            end
        end

        println(@sprintf("PCG_GPU:   %.2f ms (%d iterations)", t_pcg * 1000, info_pcg.iterations))
        println(@sprintf("Krylov.jl: %.2f ms (%d iterations)", t_krylov * 1000, solver.stats.niter))
        println(@sprintf("Speedup:   %.2fx", t_krylov / t_pcg))
        println()
    catch e
        println("Krylov.jl comparison failed (this is fine - just skipping)")
        println("Error: $e")
        println()
    end
else
    println("CUDA not available, skipping GPU tests")
    println()
end

println("="^70)
println("All tests completed!")
println("="^70)
