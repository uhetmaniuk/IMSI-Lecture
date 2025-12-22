#!/usr/bin/env julia
#
# @file PCG_GPU.jl
# @brief GPU-compatible Preconditioned Conjugate Gradient (PCG) solver
#
# Supports both CPU (SparseMatrixCSC, Vector) and GPU (CuSparseMatrixCSC, CuVector).
# Optimized for block-diagonal matrices with optional custom SpMV kernel.
#

module PCG_GPU

using LinearAlgebra
using SparseArrays

# Optional CUDA support (loaded if available)
if !isdefined(Main, :CUDA)
    try
        using CUDA
        using CUDA.CUSPARSE
        const HAS_CUDA = CUDA.functional()
    catch
        const HAS_CUDA = false
    end
else
    using CUDA
    using CUDA.CUSPARSE
    const HAS_CUDA = CUDA.functional()
end

export pcg_solve, pcg_solve!, PCGWorkspace
export extract_diagonal

"""
    PCGWorkspace{T, VecType}

Pre-allocated workspace for in-place PCG solver.
Works on both CPU (Vector) and GPU (CuVector).
"""
struct PCGWorkspace{T<:AbstractFloat, VecType<:AbstractVector{T}}
    r::VecType    # Residual
    z::VecType    # Preconditioned residual
    p::VecType    # Search direction
    Ap::VecType   # A times p
    diag::VecType # Diagonal for Jacobi preconditioner
end

"""
    PCGWorkspace(A, b)

Create workspace on the same device as the input arrays.
"""
function PCGWorkspace(A::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = size(A, 1)
    VecType = typeof(b)

    r = similar(b)
    z = similar(b)
    p = similar(b)
    Ap = similar(b)
    diag = similar(b)

    return PCGWorkspace{T, VecType}(r, z, p, Ap, diag)
end

"""
    extract_diagonal(A::SparseMatrixCSC{T}) -> Vector{T}

Extract diagonal from CPU sparse matrix.
"""
function extract_diagonal(A::SparseMatrixCSC{T}) where T
    n = size(A, 1)
    diag = zeros(T, n)

    for i = 1:n
        for k = A.colptr[i]:A.colptr[i+1]-1
            if A.rowval[k] == i
                diag[i] = A.nzval[k]
                break
            end
        end
    end

    return diag
end

# GPU diagonal extraction kernel
if HAS_CUDA
    """
    CUDA kernel to extract diagonal from CSC sparse matrix.
    Each thread handles one diagonal entry.
    """
    function extract_diagonal_kernel!(diag, colptr, rowval, nzval, n)
        i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        if i <= n
            # Search for diagonal entry in column i
            for k in colptr[i]:colptr[i+1]-1
                if rowval[k] == i
                    diag[i] = nzval[k]
                    break
                end
            end
        end

        return nothing
    end

    """
        extract_diagonal(A::CuSparseMatrixCSC{T}) -> CuVector{T}

    Extract diagonal from GPU sparse matrix using CUDA kernel.
    """
    function extract_diagonal(A::CUSPARSE.CuSparseMatrixCSC{T}) where T
        n = size(A, 1)
        diag = CUDA.zeros(T, n)

        # Launch kernel
        threads = 256
        blocks = cld(n, threads)
        @cuda threads=threads blocks=blocks extract_diagonal_kernel!(
            diag, A.colPtr, A.rowVal, A.nzVal, n
        )
        CUDA.synchronize()

        return diag
    end
end

"""
    pcg_solve(A, b; x0=nothing, tol=1e-10, maxiter=1000, verbose=false)

Preconditioned Conjugate Gradient solver with Jacobi preconditioning.
Works on both CPU and GPU.

# Arguments
- `A`: Coefficient matrix (SparseMatrixCSC or CuSparseMatrixCSC, must be SPD)
- `b`: Right-hand side vector (Vector or CuVector)
- `x0`: Initial guess (default: zero vector). **Guaranteed to be used!**
- `tol`: Convergence tolerance for ||r||_2
- `maxiter`: Maximum iterations
- `verbose`: Print iteration info

# Returns
- `x`: Solution vector (same type as b)
- `info`: NamedTuple with (iterations, converged, residual_norm)

# Example
```julia
# CPU
A_cpu = sparse(A)
b_cpu = rand(n)
x, info = pcg_solve(A_cpu, b_cpu; x0=x_initial, tol=1e-12)

# GPU
A_gpu = CuSparseMatrixCSC(A)
b_gpu = CuArray(b)
x_gpu, info = pcg_solve(A_gpu, b_gpu; x0=CuArray(x_initial), tol=1e-12)
```
"""
function pcg_solve(A::AbstractMatrix{T}, b::AbstractVector{T};
                   x0::Union{Nothing, AbstractVector{T}}=nothing,
                   tol::T=T(1e-10), maxiter::Int=1000,
                   verbose::Bool=false) where T<:AbstractFloat

    n = length(b)
    @assert size(A) == (n, n) "Matrix A must be square and match size of b"

    # Create workspace
    workspace = PCGWorkspace(A, b)

    # Extract diagonal for Jacobi preconditioner
    extract_diagonal!(workspace.diag, A)

    # Initialize solution
    x = similar(b)
    if x0 !== nothing
        copyto!(x, x0)  # Use initial guess
    else
        fill!(x, zero(T))  # Start from zero
    end

    # Call in-place solver
    info = pcg_solve!(x, workspace, A, b; tol=tol, maxiter=maxiter, verbose=verbose)

    return x, info
end

"""
    extract_diagonal!(diag, A)

Extract diagonal into pre-allocated array (dispatches to CPU or GPU version).
"""
function extract_diagonal!(diag::AbstractVector{T}, A::AbstractMatrix{T}) where T
    d = extract_diagonal(A)
    copyto!(diag, d)
    return diag
end

"""
    pcg_solve!(x, workspace, A, b; tol=1e-10, maxiter=1000, verbose=false)

In-place PCG solver. Modifies x in-place.
Uses pre-allocated workspace to avoid allocations.

# Note
- `x` serves as both initial guess and output
- Workspace diagonal must be pre-filled with extract_diagonal!(workspace.diag, A)
- Zero allocations during solve
"""
function pcg_solve!(x::AbstractVector{T}, workspace::PCGWorkspace{T},
                   A::AbstractMatrix{T}, b::AbstractVector{T};
                   tol::T=T(1e-10), maxiter::Int=1000,
                   verbose::Bool=false) where T<:AbstractFloat

    # Unpack workspace
    r = workspace.r
    z = workspace.z
    p = workspace.p
    Ap = workspace.Ap
    diag = workspace.diag

    n = length(b)

    # Compute initial residual: r = b - A*x
    mul!(r, A, x)           # r = A*x
    r .= b .- r             # r = b - r

    # Apply Jacobi preconditioner: z = r ./ diag
    # Handle zero diagonal entries
    @. z = ifelse(diag != 0, r / diag, r)

    # Initialize search direction: p = z
    copyto!(p, z)

    # Compute initial (r, z)
    rz_old = dot(r, z)

    # Track initial residual
    initial_residual = norm(b)
    if verbose
        println("PCG: initial ||b|| = ", initial_residual)
    end

    # PCG iteration
    for iter = 1:maxiter
        # Compute Ap = A * p
        mul!(Ap, A, p)

        # Compute (Ap, p)
        pAp = dot(Ap, p)

        # Check for breakdown (matrix not positive definite)
        if pAp <= 0
            if verbose
                @warn "PCG breakdown: matrix not positive definite at iteration $iter"
            end
            return (iterations=-1, converged=false, residual_norm=T(Inf))
        end

        # Compute step length: alpha = (r, z) / (Ap, p)
        alpha = rz_old / pAp

        # Update solution and residual
        # x += alpha * p
        # r -= alpha * Ap
        x .+= alpha .* p
        r .-= alpha .* Ap

        # Check convergence
        residual_norm = norm(r)

        if verbose && (iter % 10 == 0 || residual_norm < tol)
            relative = residual_norm / max(initial_residual, eps(T))
            @printf("  PCG iter %4d: ||r|| = %.2e (rel = %.2e)\n",
                    iter, residual_norm, relative)
        end

        if residual_norm < tol
            if verbose
                println("PCG converged in $iter iterations")
            end
            return (iterations=iter, converged=true, residual_norm=residual_norm)
        end

        # Apply Jacobi preconditioner: z = r ./ diag
        @. z = ifelse(diag != 0, r / diag, r)

        # Compute new (r, z)
        rz_new = dot(r, z)

        # Compute beta
        beta = rz_new / rz_old

        # Update search direction: p = z + beta * p
        p .= z .+ beta .* p

        # Update (r, z) for next iteration
        rz_old = rz_new
    end

    # Did not converge
    residual_norm = norm(r)
    if verbose
        @warn "PCG did not converge in $maxiter iterations (residual = $residual_norm)"
    end

    return (iterations=-1, converged=false, residual_norm=residual_norm)
end

# Helper for formatted printing (add Printf if not already loaded)
using Printf

end # module PCG_GPU
