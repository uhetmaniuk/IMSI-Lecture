using CUDA
using BenchmarkTools
using SparseArrays
using LinearAlgebra

# =============================================================================
# Mesh: nx×ny uniform rectangles on [0,1]²  (same as previous files)
#   Node numbering (0-based i,j): node(i,j) = j*(nx+1) + i + 1
#   Q1 local node order: 1=SW, 2=SE, 3=NE, 4=NW
# =============================================================================
function make_mesh(nx, ny)
    nnodes = (nx+1)*(ny+1)
    nelems = nx * ny
    coords = Matrix{Float64}(undef, 2, nnodes)
    conn   = Matrix{Int32}(undef, 4, nelems)    # Int32: halves bandwidth on GPU
    colors = Vector{Int}(undef, nelems)
    for j in 0:ny, i in 0:nx
        n = j*(nx+1) + i + 1
        coords[1, n] = i / nx
        coords[2, n] = j / ny
    end
    for j in 0:ny-1, i in 0:nx-1
        sw = j*(nx+1) + i + 1
        e  = j*nx + i + 1
        conn[1, e] = sw
        conn[2, e] = sw + 1
        conn[3, e] = sw + (nx+1) + 1
        conn[4, e] = sw + (nx+1)
        colors[e]  = 1 + (i & 1) + 2*(j & 1)
    end
    nnodes, coords, conn, colors
end

# =============================================================================
# Reference: assemble sparse K on CPU, then use mul! for SpMV.
#   Used only to verify the matrix-free results.
# =============================================================================
function assemble_K_cpu(coords, conn, nnodes)
    nelems = size(conn, 2)
    I = Vector{Int}(undef, 16*nelems)
    J = Vector{Int}(undef, 16*nelems)
    V = Vector{Float64}(undef, 16*nelems)
    g = 1.0 / sqrt(3.0)
    k = 1
    for e in 1:nelems
        n = (Int(conn[1,e]), Int(conn[2,e]), Int(conn[3,e]), Int(conn[4,e]))
        x = (coords[1,n[1]], coords[1,n[2]], coords[1,n[3]], coords[1,n[4]])
        y = (coords[2,n[1]], coords[2,n[2]], coords[2,n[3]], coords[2,n[4]])
        Ke = zeros(4, 4)
        for η in (-g, g), ξ in (-g, g)
            dNdξ = (-(1-η), (1-η), (1+η), -(1+η)) .* 0.25
            dNdη = (-(1-ξ), -(1+ξ), (1+ξ), (1-ξ)) .* 0.25
            J11 = sum(x .* dNdξ); J12 = sum(x .* dNdη)
            J21 = sum(y .* dNdξ); J22 = sum(y .* dNdη)
            detJ = J11*J22 - J12*J21; iJ = 1.0/detJ
            dNdx = ( J22.*dNdξ .- J21.*dNdη) .* iJ
            dNdy = (-J12.*dNdξ .+ J11.*dNdη) .* iJ
            Ke .+= detJ .* (dNdx .* dNdx' .+ dNdy .* dNdy')
        end
        for b in 1:4, a in 1:4
            I[k] = n[a]; J[k] = n[b]; V[k] = Ke[a,b]; k += 1
        end
    end
    sparse(I, J, V, nnodes, nnodes)
end

# =============================================================================
# CPU matrix-free SpMV  (single-threaded baseline)
#
#   For each element: gather u, compute Ke·u_loc on the fly, scatter.
#
#   Key insight: instead of forming Ke (4×4) and then multiplying,
#   we compute B·u_loc first (a length-2 vector at each Gauss point)
#   and then B'*(B·u_loc).  This reduces the inner work from
#   O(4²) = 16 FMAs for Ke storage + 4×4 multiply
#   to O(4+4) = 8 FMAs per Gauss point — half the FLOPs, zero extra storage.
# =============================================================================
function matvec_cpu!(v, u, coords, conn)
    fill!(v, 0.0)
    g = 1.0 / sqrt(3.0)
    for e in axes(conn, 2)
        n1, n2, n3, n4 = conn[1,e], conn[2,e], conn[3,e], conn[4,e]
        x1, y1 = coords[1,n1], coords[2,n1]
        x2, y2 = coords[1,n2], coords[2,n2]
        x3, y3 = coords[1,n3], coords[2,n3]
        x4, y4 = coords[1,n4], coords[2,n4]
        u1, u2, u3, u4 = u[n1], u[n2], u[n3], u[n4]
        v1 = v2 = v3 = v4 = 0.0
        for η in (-g, g), ξ in (-g, g)
            dN1dξ = -(1-η)*0.25; dN2dξ =  (1-η)*0.25
            dN3dξ =  (1+η)*0.25; dN4dξ = -(1+η)*0.25
            dN1dη = -(1-ξ)*0.25; dN2dη = -(1+ξ)*0.25
            dN3dη =  (1+ξ)*0.25; dN4dη =  (1-ξ)*0.25
            J11 = x1*dN1dξ + x2*dN2dξ + x3*dN3dξ + x4*dN4dξ
            J12 = x1*dN1dη + x2*dN2dη + x3*dN3dη + x4*dN4dη
            J21 = y1*dN1dξ + y2*dN2dξ + y3*dN3dξ + y4*dN4dξ
            J22 = y1*dN1dη + y2*dN2dη + y3*dN3dη + y4*dN4dη
            detJ = J11*J22 - J12*J21;  iJ = 1.0 / detJ
            dN1dx = ( J22*dN1dξ - J21*dN1dη)*iJ; dN1dy = (-J12*dN1dξ + J11*dN1dη)*iJ
            dN2dx = ( J22*dN2dξ - J21*dN2dη)*iJ; dN2dy = (-J12*dN2dξ + J11*dN2dη)*iJ
            dN3dx = ( J22*dN3dξ - J21*dN3dη)*iJ; dN3dy = (-J12*dN3dξ + J11*dN3dη)*iJ
            dN4dx = ( J22*dN4dξ - J21*dN4dη)*iJ; dN4dy = (-J12*dN4dξ + J11*dN4dη)*iJ
            # B·u: project u onto gradient space at this Gauss point
            Bu_x = dN1dx*u1 + dN2dx*u2 + dN3dx*u3 + dN4dx*u4
            Bu_y = dN1dy*u1 + dN2dy*u2 + dN3dy*u3 + dN4dy*u4
            # B'*(B·u): scatter back, weighted by detJ (Gauss weight = 1)
            v1 += detJ * (dN1dx*Bu_x + dN1dy*Bu_y)
            v2 += detJ * (dN2dx*Bu_x + dN2dy*Bu_y)
            v3 += detJ * (dN3dx*Bu_x + dN3dy*Bu_y)
            v4 += detJ * (dN4dx*Bu_x + dN4dy*Bu_y)
        end
        v[n1] += v1;  v[n2] += v2;  v[n3] += v3;  v[n4] += v4
    end
    v
end

# =============================================================================
# GPU kernel — variant 1: one thread per element, atomic scatter
#
#   Each thread independently:
#     1. reads 8 coords + 4 u values  (gather)
#     2. runs the 2×2 Gauss loop      (compute Ke·u_loc on the fly)
#     3. writes 4 results with atomicAdd (scatter)
#
#   atomicAdd is needed because multiple elements share nodes.
#   On modern GPUs (Volta+) Float64 atomics are hardware-accelerated.
#   The compute phase is long enough to hide most of the atomic latency.
# =============================================================================
function kernel_matvec_atomic!(v, u, coords, conn, nelems)
    e = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    e > nelems && return nothing

    n1 = conn[1,e];  n2 = conn[2,e];  n3 = conn[3,e];  n4 = conn[4,e]
    x1 = coords[1,n1];  y1 = coords[2,n1]
    x2 = coords[1,n2];  y2 = coords[2,n2]
    x3 = coords[1,n3];  y3 = coords[2,n3]
    x4 = coords[1,n4];  y4 = coords[2,n4]
    u1 = u[n1];  u2 = u[n2];  u3 = u[n3];  u4 = u[n4]

    g = 1.0 / sqrt(3.0)
    v1 = v2 = v3 = v4 = 0.0

    for η in (-g, g), ξ in (-g, g)
        dN1dξ = -(1-η)*0.25; dN2dξ =  (1-η)*0.25
        dN3dξ =  (1+η)*0.25; dN4dξ = -(1+η)*0.25
        dN1dη = -(1-ξ)*0.25; dN2dη = -(1+ξ)*0.25
        dN3dη =  (1+ξ)*0.25; dN4dη =  (1-ξ)*0.25
        J11 = x1*dN1dξ + x2*dN2dξ + x3*dN3dξ + x4*dN4dξ
        J12 = x1*dN1dη + x2*dN2dη + x3*dN3dη + x4*dN4dη
        J21 = y1*dN1dξ + y2*dN2dξ + y3*dN3dξ + y4*dN4dξ
        J22 = y1*dN1dη + y2*dN2dη + y3*dN3dη + y4*dN4dη
        detJ = J11*J22 - J12*J21;  iJ = 1.0 / detJ
        dN1dx = ( J22*dN1dξ - J21*dN1dη)*iJ; dN1dy = (-J12*dN1dξ + J11*dN1dη)*iJ
        dN2dx = ( J22*dN2dξ - J21*dN2dη)*iJ; dN2dy = (-J12*dN2dξ + J11*dN2dη)*iJ
        dN3dx = ( J22*dN3dξ - J21*dN3dη)*iJ; dN3dy = (-J12*dN3dξ + J11*dN3dη)*iJ
        dN4dx = ( J22*dN4dξ - J21*dN4dη)*iJ; dN4dy = (-J12*dN4dξ + J11*dN4dη)*iJ
        Bu_x = dN1dx*u1 + dN2dx*u2 + dN3dx*u3 + dN4dx*u4
        Bu_y = dN1dy*u1 + dN2dy*u2 + dN3dy*u3 + dN4dy*u4
        v1 += detJ * (dN1dx*Bu_x + dN1dy*Bu_y)
        v2 += detJ * (dN2dx*Bu_x + dN2dy*Bu_y)
        v3 += detJ * (dN3dx*Bu_x + dN3dy*Bu_y)
        v4 += detJ * (dN4dx*Bu_x + dN4dy*Bu_y)
    end

    # atomicAdd: thread-safe scatter without locks.
    # Safe to use Float64 on Volta (sm_70) and later.
    CUDA.@atomic v[n1] += v1
    CUDA.@atomic v[n2] += v2
    CUDA.@atomic v[n3] += v3
    CUDA.@atomic v[n4] += v4

    return nothing
end

# Wrapper: allocate nothing, just configure and launch the kernel.
function matvec_gpu_atomic!(v_d, u_d, coords_d, conn_d, nelems)
    fill!(v_d, 0.0)
    nthreads = 256
    nblocks  = cld(nelems, nthreads)
    @cuda threads=nthreads blocks=nblocks kernel_matvec_atomic!(
        v_d, u_d, coords_d, conn_d, nelems)
    v_d
end

# =============================================================================
# GPU kernel — variant 2: coloring, no atomics
#
#   Process one color group at a time.  Within a group, no two elements
#   share a node, so threads can write directly to v without atomicAdd.
#   The host launches 4 kernels sequentially (one per color); each kernel
#   processes roughly nelems/4 elements.
#
#   Trade-off vs. atomic variant:
#     + No atomic contention → better throughput on older GPUs (pre-Volta)
#     + Writes coalesce better within a warp if elements in a group are
#       laid out contiguously (which they are after grouping by color)
#     - 4 kernel launches instead of 1; launch overhead ~5–10 µs each
#     - Requires storing color_group index arrays on the device
# =============================================================================
function kernel_matvec_color!(v, u, coords, conn, group, ngroup)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    idx > ngroup && return nothing
    e = group[idx]

    n1 = conn[1,e];  n2 = conn[2,e];  n3 = conn[3,e];  n4 = conn[4,e]
    x1 = coords[1,n1];  y1 = coords[2,n1]
    x2 = coords[1,n2];  y2 = coords[2,n2]
    x3 = coords[1,n3];  y3 = coords[2,n3]
    x4 = coords[1,n4];  y4 = coords[2,n4]
    u1 = u[n1];  u2 = u[n2];  u3 = u[n3];  u4 = u[n4]

    g = 1.0 / sqrt(3.0)
    v1 = v2 = v3 = v4 = 0.0

    for η in (-g, g), ξ in (-g, g)
        dN1dξ = -(1-η)*0.25; dN2dξ =  (1-η)*0.25
        dN3dξ =  (1+η)*0.25; dN4dξ = -(1+η)*0.25
        dN1dη = -(1-ξ)*0.25; dN2dη = -(1+ξ)*0.25
        dN3dη =  (1+ξ)*0.25; dN4dη =  (1-ξ)*0.25
        J11 = x1*dN1dξ + x2*dN2dξ + x3*dN3dξ + x4*dN4dξ
        J12 = x1*dN1dη + x2*dN2dη + x3*dN3dη + x4*dN4dη
        J21 = y1*dN1dξ + y2*dN2dξ + y3*dN3dξ + y4*dN4dξ
        J22 = y1*dN1dη + y2*dN2dη + y3*dN3dη + y4*dN4dη
        detJ = J11*J22 - J12*J21;  iJ = 1.0 / detJ
        dN1dx = ( J22*dN1dξ - J21*dN1dη)*iJ; dN1dy = (-J12*dN1dξ + J11*dN1dη)*iJ
        dN2dx = ( J22*dN2dξ - J21*dN2dη)*iJ; dN2dy = (-J12*dN2dξ + J11*dN2dη)*iJ
        dN3dx = ( J22*dN3dξ - J21*dN3dη)*iJ; dN3dy = (-J12*dN3dξ + J11*dN3dη)*iJ
        dN4dx = ( J22*dN4dξ - J21*dN4dη)*iJ; dN4dy = (-J12*dN4dξ + J11*dN4dη)*iJ
        Bu_x = dN1dx*u1 + dN2dx*u2 + dN3dx*u3 + dN4dx*u4
        Bu_y = dN1dy*u1 + dN2dy*u2 + dN3dy*u3 + dN4dy*u4
        v1 += detJ * (dN1dx*Bu_x + dN1dy*Bu_y)
        v2 += detJ * (dN2dx*Bu_x + dN2dy*Bu_y)
        v3 += detJ * (dN3dx*Bu_x + dN3dy*Bu_y)
        v4 += detJ * (dN4dx*Bu_x + dN4dy*Bu_y)
    end

    # Direct write — safe because same-color elements do not share nodes.
    v[n1] += v1;  v[n2] += v2;  v[n3] += v3;  v[n4] += v4

    return nothing
end

function matvec_gpu_color!(v_d, u_d, coords_d, conn_d, color_groups_d)
    fill!(v_d, 0.0)
    nthreads = 256
    for group_d in color_groups_d
        ngroup = length(group_d)
        nblocks = cld(ngroup, nthreads)
        @cuda threads=nthreads blocks=nblocks kernel_matvec_color!(
            v_d, u_d, coords_d, conn_d, group_d, ngroup)
        # Implicit device-side barrier: the next @cuda call cannot begin
        # until the current kernel finishes (within the same CUDA stream).
    end
    v_d
end

# =============================================================================
# Driver
# =============================================================================
nx, ny = 1024, 1024
nnodes, coords, conn, colors = make_mesh(nx, ny)
nelems = nx * ny
color_groups = [findall(==(c), colors) for c in 1:4]

println("Mesh : $(nx)×$(ny), $nnodes nodes, $nelems elements")
println("GPU  : ", CUDA.name(CUDA.device()))
println()

# Random test vector
u_cpu = randn(Float64, nnodes)

# --- Reference: assembled sparse K on CPU ---
print("Assembling reference K … ")
@time K_ref = assemble_K_cpu(coords, conn, nnodes)
v_ref = K_ref * u_cpu

# --- CPU matrix-free ---
v_cpu = zeros(Float64, nnodes)
matvec_cpu!(v_cpu, u_cpu, coords, conn)   # warm-up
println("Max error (CPU matrix-free) : ", maximum(abs.(v_cpu .- v_ref)))

# --- Upload to GPU ---
coords_d       = CuArray(coords)
conn_d         = CuArray(conn)
u_d            = CuArray(u_cpu)
v_d            = CUDA.zeros(Float64, nnodes)
color_groups_d = [CuArray(Int32.(g)) for g in color_groups]

# --- GPU atomic ---
matvec_gpu_atomic!(v_d, u_d, coords_d, conn_d, nelems)   # warm-up
v_gpu_atomic = Array(v_d)
println("Max error (GPU atomic)      : ", maximum(abs.(v_gpu_atomic .- v_ref)))

# --- GPU coloring ---
matvec_gpu_color!(v_d, u_d, coords_d, conn_d, color_groups_d)   # warm-up
v_gpu_color = Array(v_d)
println("Max error (GPU coloring)    : ", maximum(abs.(v_gpu_color .- v_ref)))
println()

# --- Benchmarks ---
# CUDA.@sync ensures the GPU work is complete before the timer stops.
print("CPU sparse K*u (ref)  : ")
@btime mul!($v_cpu, $K_ref, $u_cpu)

print("CPU matrix-free       : ")
@btime matvec_cpu!($v_cpu, $u_cpu, $coords, $conn)

print("GPU atomic            : ")
@btime CUDA.@sync matvec_gpu_atomic!($v_d, $u_d, $coords_d, $conn_d, $nelems)

print("GPU coloring          : ")
@btime CUDA.@sync matvec_gpu_color!($v_d, $u_d, $coords_d, $conn_d, $color_groups_d)

