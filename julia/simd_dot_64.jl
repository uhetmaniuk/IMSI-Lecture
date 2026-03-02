using BenchmarkTools
using Printf
using LinearAlgebra

function dot_naive(a::Vector{Float64}, b::Vector{Float64})
    s = Float64(0)
    for i in eachindex(a)
        @inbounds s += Float64(a[i]) * Float64(b[i])
    end
    return Float64(s)
end

function dot_simd(a::Vector{Float64}, b::Vector{Float64})
    s = Float64(0)
    @simd for i in eachindex(a)
        s += a[i] * b[i]
    end
    return Float64(s)
end

p = 10 .^ (2:0.1:6)

for n in ceil.(Int, p)
a = rand(n)
b = rand(n)

snaive = dot_naive(a, b)
ssimd = dot_simd(a, b)
sref = dot(a, b)

println(" Test ", @sprintf("%d %d %.6e %.6e", BLAS.get_num_threads(), n, (snaive - sref)/sref, (ssimd - sref)/sref))

println("Naive:")
@btime dot_naive($a, $b)

println("With @simd:")
@btime dot_simd($a, $b)

println("Built-in dot (for reference):")
using LinearAlgebra
@btime dot($a, $b)
end
