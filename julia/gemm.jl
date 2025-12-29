using LinearAlgebra
using Printf

"""
    gemm!(C, A, B)

Compute C += A * B using threaded outer product algorithm.
Parallelizes over columns of B.
"""
function gemm!(C, A, B)
    A_rows, A_cols = size(A)
    B_rows, B_cols = size(B)

    @assert A_cols == B_rows "Inner dimensions must match"
    @assert size(C) == (A_rows, B_cols) "Output size mismatch"

    # With @inbounds (unsafe) - skips the check on array indices
    # Bounds checking has overhead (small but adds up in tight loops)
    Threads.@threads for j in 1:B_cols
        for l in 1:A_cols
            @inbounds temp = B[l, j]
            for i in 1:A_rows
                @inbounds C[i, j] += temp * A[i, l]
            end
        end
    end

    return C
end

"""
    benchmark_gemm(n, nthreads, nruns=3)

Benchmark matrix multiplication for n×n matrices with specified thread count.
Returns average time in seconds.
"""
function benchmark_gemm(n::Int, nthreads::Int, nruns::Int=3)
    A = rand(n, n)
    B = rand(n, n)

    # Set thread count (note: this affects subsequent operations)
    # Julia doesn't allow dynamic thread count changes, so this is informational

    # Warmup
    C_test = zeros(n, n)
    gemm!(C_test, A, B)

    # Benchmark custom implementation
    times = Float64[]
    for _ in 1:nruns
        C = zeros(n, n)
        t = @elapsed gemm!(C, A, B)
        push!(times, t)
    end

    return minimum(times)
end

"""
    run_gemm_tests()

Run GEMM benchmark for various matrix sizes with current thread count.
"""
function run_gemm_tests()
    nthreads = Threads.nthreads()

    # Test sizes
    sizes = [64, 128, 256, 512, 1024, 2048]

    for n in sizes
        time = benchmark_gemm(n, nthreads, 3)
        @printf("%d,%d,%.6f\n", n, nthreads, time)
        flush(stdout)
    end
end

# Run tests if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_gemm_tests()
end
