using CUDA
using BenchmarkTools

# ============================================================================
# Block-level reduction using cooperative groups and shared memory
# ============================================================================

"""
    block_reduce_sum(g::CG.ThreadBlock, shmem, val)

Perform a sum reduction across all threads in a block.
Returns the reduced value on thread 1 of the block.

# Arguments
- `g`: Thread block group from CG.this_thread_block()
- `shmem`: Shared memory array (CuStaticSharedArray or CuDynamicSharedArray)
- `val`: This thread's value to contribute to the reduction

# Returns
The sum of all values across the block (valid only on thread 1)
"""
function block_reduce_sum(g, shmem, val)
    tid = threadIdx().x
    
    # Store this thread's value in shared memory
    shmem[tid] = val
    CG.sync(g)
    
    # Parallel tree reduction in shared memory
    s = blockDim().x ÷ 2
    while s > 0
        if tid <= s
            shmem[tid] += shmem[tid + s]
        end
        CG.sync(g)
        s ÷= 2
    end
    
    # Result is in shmem[1]
    return shmem[1]
end

"""
    block_reduce_general(g::CG.ThreadBlock, shmem, val, op)

Perform a reduction across all threads in a block with a custom binary operator.

# Arguments
- `g`: Thread block group
- `shmem`: Shared memory array
- `val`: This thread's value
- `op`: Binary reduction operator (e.g., +, *, min, max)
"""
function block_reduce_general(g, shmem, val, op)
    tid = threadIdx().x
    
    shmem[tid] = val
    CG.sync(g)
    
    s = blockDim().x ÷ 2
    while s > 0
        if tid <= s
            shmem[tid] = op(shmem[tid], shmem[tid + s])
        end
        CG.sync(g)
        s ÷= 2
    end
    
    return shmem[1]
end

# ============================================================================
# Example 1: Simple sum reduction kernel
# ============================================================================

function sum_reduce_kernel(input::CuDeviceArray{T}, output::CuDeviceArray{T}) where T
    # Get thread block group
    g = CG.this_thread_block()
    
    # Allocate shared memory (must match block size)
    shmem = CuStaticSharedArray(T, 256)
    
    # Each thread loads one element (or zero if out of bounds)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    val = idx <= length(input) ? input[idx] : zero(T)
    
    # Perform block-level reduction
    block_sum = block_reduce_sum(g, shmem, val)
    
    # First thread writes result
    if threadIdx().x == 1
        output[blockIdx().x] = block_sum
    end
    
    return nothing
end

# ============================================================================
# Example 2: Reduction with multiple elements per thread (grid-stride)
# ============================================================================

function sum_reduce_multielement_kernel(input::CuDeviceArray{T}, 
                                        output::CuDeviceArray{T}) where T
    g = CG.this_thread_block()
    shmem = CuStaticSharedArray(T, 256)
    
    # Global thread index
    tid = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    
    # Each thread accumulates multiple elements
    thread_sum = zero(T)
    idx = tid
    while idx <= length(input)
        thread_sum += input[idx]
        idx += stride
    end
    
    # Reduce across the block
    block_sum = block_reduce_sum(g, shmem, thread_sum)
    
    if threadIdx().x == 1
        CUDA.@atomic output[1] += block_sum
    end
    
    return nothing
end

# ============================================================================
# Example 3: Max reduction
# ============================================================================

function max_reduce_kernel(input::CuDeviceArray{T}, output::CuDeviceArray{T}) where T
    g = CG.this_thread_block()
    shmem = CuStaticSharedArray(T, 256)
    
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    val = idx <= length(input) ? input[idx] : typemin(T)
    
    block_max = block_reduce_general(g, shmem, val, max)
    
    if threadIdx().x == 1
        CUDA.@atomic output[1] = max(output[1], block_max)
    end
    
    return nothing
end

# ============================================================================
# Example 4: Dot product (reduction of products)
# ============================================================================

function dot_product_kernel(a::CuDeviceArray{T}, b::CuDeviceArray{T}, 
                           output::CuDeviceArray{T}) where T
    g = CG.this_thread_block()
    shmem = CuStaticSharedArray(T, 256)
    
    # Global thread index with grid-stride loop
    tid = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    
    # Compute partial dot product
    thread_dot = zero(T)
    idx = tid
    while idx <= length(a)
        thread_dot += a[idx] * b[idx]
        idx += stride
    end
    
    # Reduce across block
    block_dot = block_reduce_sum(g, shmem, thread_dot)
    
    if threadIdx().x == 1
        CUDA.@atomic output[1] += block_dot
    end
    
    return nothing
end

# ============================================================================
# Example 5: Using dynamic shared memory for flexible block sizes
# ============================================================================

function sum_reduce_dynamic_shmem_kernel(input::CuDeviceArray{T}, 
                                         output::CuDeviceArray{T}) where T
    g = CG.this_thread_block()
    
    # Use dynamic shared memory (size specified at launch time)
    shmem = CuDynamicSharedArray(T, blockDim().x)
    
    tid = threadIdx().x
    idx = (blockIdx().x - 1) * blockDim().x + tid
    val = idx <= length(input) ? input[idx] : zero(T)
    
    # Reduction
    shmem[tid] = val
    CG.sync(g)
    
    s = blockDim().x ÷ 2
    while s > 0
        if tid <= s
            shmem[tid] += shmem[tid + s]
        end
        CG.sync(g)
        s ÷= 2
    end
    
    if tid == 1
        output[blockIdx().x] = shmem[1]
    end
    
    return nothing
end

# ============================================================================
# Helper functions and tests
# ============================================================================

function test_sum_reduction()
    println("\n=== Testing Sum Reduction ===")
    
    n = 1024 * 1024
    data = CUDA.rand(Float32, n)
    
    # CPU reference
    cpu_sum = sum(Array(data))
    println("CPU sum: $cpu_sum")
    
    # GPU reduction
    threads = 256
    blocks = cld(n, threads)
    
    # Method 1: Two-stage reduction
    partial_sums = CUDA.zeros(Float32, blocks)
    @cuda threads=threads blocks=blocks sum_reduce_kernel(data, partial_sums)
    gpu_sum = sum(Array(partial_sums))
    println("GPU sum (two-stage): $gpu_sum")
    println("Relative error: $(abs(gpu_sum - cpu_sum) / cpu_sum)")
    
    # Method 2: Single atomic reduction
    result = CUDA.zeros(Float32, 1)
    @cuda threads=threads blocks=min(blocks, 1024) sum_reduce_multielement_kernel(data, result)
    gpu_sum2 = Array(result)[1]
    println("GPU sum (atomic): $gpu_sum2")
    println("Relative error: $(abs(gpu_sum2 - cpu_sum) / cpu_sum)")
end

function test_max_reduction()
    println("\n=== Testing Max Reduction ===")
    
    n = 1024 * 1024
    data = CUDA.rand(Float32, n)
    
    cpu_max = maximum(Array(data))
    println("CPU max: $cpu_max")
    
    result = CUDA.fill(typemin(Float32), 1)
    threads = 256
    blocks = min(cld(n, threads), 1024)
    
    @cuda threads=threads blocks=blocks max_reduce_kernel(data, result)
    gpu_max = Array(result)[1]
    println("GPU max: $gpu_max")
    println("Difference: $(abs(gpu_max - cpu_max))")
end

function test_dot_product()
    println("\n=== Testing Dot Product ===")
    
    n = 1024 * 1024
    a = CUDA.rand(Float32, n)
    b = CUDA.rand(Float32, n)
    
    cpu_dot = sum(Array(a) .* Array(b))
    println("CPU dot: $cpu_dot")
    
    result = CUDA.zeros(Float32, 1)
    threads = 256
    blocks = min(cld(n, threads), 1024)
    
    @cuda threads=threads blocks=blocks dot_product_kernel(a, b, result)
    gpu_dot = Array(result)[1]
    println("GPU dot: $gpu_dot")
    println("Relative error: $(abs(gpu_dot - cpu_dot) / cpu_dot)")
end

function test_dynamic_shmem()
    println("\n=== Testing Dynamic Shared Memory ===")
    
    n = 2048
    data = CUDA.rand(Float32, n)
    
    cpu_sum = sum(Array(data))
    
    # Test with different block sizes
    for threads in [64, 128, 256, 512]
        blocks = cld(n, threads)
        partial_sums = CUDA.zeros(Float32, blocks)
        
        # Launch with dynamic shared memory
        @cuda threads=threads blocks=blocks shmem=threads*sizeof(Float32) sum_reduce_dynamic_shmem_kernel(data, partial_sums)
        
        gpu_sum = sum(Array(partial_sums))
        err = abs(gpu_sum - cpu_sum) / cpu_sum
        println("Block size $threads: sum=$gpu_sum, rel_error=$err")
    end
end

function benchmark_reductions()
    println("\n=== Benchmarking ===")
    
    n = 8 * 1024 * 1024
    data = CUDA.rand(Float32, n)
    
    threads = 256
    blocks = min(cld(n, threads), 1024)
    
    # Warm up
    result = CUDA.zeros(Float32, 1)
    @cuda threads=threads blocks=blocks sum_reduce_multielement_kernel(data, result)
    CUDA.synchronize()
    
    # Benchmark
    println("Array size: $n elements")
    t = CUDA.@elapsed begin
        result = CUDA.zeros(Float32, 1)
        @cuda threads=threads blocks=blocks sum_reduce_multielement_kernel(data, result)
        CUDA.synchronize()
    end
    
    bandwidth = (n * sizeof(Float32)) / t / 1e9
    println("Time: $(t*1000) ms")
    println("Bandwidth: $(bandwidth) GB/s")
end

# ============================================================================
# Main execution
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    println("CUDA.jl Block-Level Reduction Examples")
    println("=" ^ 50)
    
    test_sum_reduction()
    test_max_reduction()
    test_dot_product()
    test_dynamic_shmem()
    benchmark_reductions()
end
