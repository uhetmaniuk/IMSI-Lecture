using Base.Threads
using BenchmarkTools
using StaticArrays
using SparseArrays
using LinearAlgebra

# =============================================================================
# Mesh: nx×ny uniform rectangles on [0,1]²
#   Node numbering (0-based i,j): node(i,j) = j*(nx+1) + i + 1
#   Q1 local node order: 1=SW, 2=SE, 3=NE, 4=NW
#   coords[2, nnodes]: column k holds (x,y) of node k
# =============================================================================
function make_mesh(nx, ny)
    nnodes = (nx+1)*(ny+1)
    nelems = nx * ny
    coords = Matrix{Float64}(undef, 2, nnodes)
    conn   = Matrix{Int}(undef, 4, nelems)
    colors = Vector{Int}(undef, nelems)

    for j in 0:ny, i in 0:nx
        n = j*(nx+1) + i + 1
        coords[1, n] = i / nx
        coords[2, n] = j / ny
    end

    for j in 0:ny-1, i in 0:nx-1
        sw = j*(nx+1) + i + 1
        e  = j*nx + i + 1
        conn[1, e] = sw                  # SW
        conn[2, e] = sw + 1              # SE
        conn[3, e] = sw + (nx+1) + 1    # NE
        conn[4, e] = sw + (nx+1)         # NW
        colors[e]  = 1 + (i & 1) + 2*(j & 1)   # 4-color checkerboard
    end

    nnodes, coords, conn, colors
end

# =============================================================================
# Compute the 4×4 element stiffness matrix for Laplace given the four
# corner coordinates of one Q1 rectangle.
#   2×2 Gauss quadrature (exact for Q1 on parallelogram elements).
#   xy: 2×4 matrix, columns = (x,y) of nodes 1..4 in local order.
# =============================================================================
function element_Ke(xy::SMatrix{2,4,Float64,8})
    g  = 1.0 / sqrt(3.0)
    Ke = zeros(MMatrix{4,4,Float64,16})

    for η in (-g, g), ξ in (-g, g)
        # Shape function derivatives in reference space
        dNdξ = SVector{4}( -(1-η),  (1-η),  (1+η), -(1+η) ) * 0.25
        dNdη = SVector{4}( -(1-ξ), -(1+ξ),  (1+ξ),  (1-ξ) ) * 0.25

        # Jacobian J[2×2] = xy[2×4] * [dNdξ dNdη][4×2]
        J    = xy * hcat(dNdξ, dNdη)   # SMatrix{2,2}
        J11, J12, J21, J22 = J[1,1], J[1,2], J[2,1], J[2,2]
        detJ = J11*J22 - J12*J21

        # Physical gradients via J⁻¹
        inv_detJ = 1.0 / detJ
        dNdx = ( J22 .* dNdξ .- J21 .* dNdη) .* inv_detJ
        dNdy = (-J12 .* dNdξ .+ J11 .* dNdη) .* inv_detJ

        Ke .+= detJ .* (dNdx * dNdx' .+ dNdy * dNdy')
    end
    SMatrix(Ke)
end

# Convenience wrapper: gather node coordinates from global arrays and call above
@inline function element_Ke(coords, nodes)
    xy = SMatrix{2,4,Float64,8}(
        coords[1,nodes[1]], coords[2,nodes[1]],
        coords[1,nodes[2]], coords[2,nodes[2]],
        coords[1,nodes[3]], coords[2,nodes[3]],
        coords[1,nodes[4]], coords[2,nodes[4]],
    )
    element_Ke(xy)
end

# =============================================================================
# Step 1 — Build sparsity pattern and the elem_ptrs index table.
#
#   elem_ptrs[16, nelems]: for each element e and each of its 16 local (a,b)
#   pairs (a=row, b=col, both in 1:4), stores the position in K.nzval where
#   K[conn[a,e], conn[b,e]] lives.
#
#   Building elem_ptrs once means assembly never has to search; it only does
#   direct indexed writes into nzval, which is safe under coloring.
# =============================================================================
function build_sparse_K(coords, conn, nnodes)
    nelems = size(conn, 2)

    # --- Collect all (row, col) pairs from element connectivities ---
    #     Each element contributes 4×4 = 16 pairs.
    npairs = 16 * nelems
    I = Vector{Int}(undef, npairs)
    J = Vector{Int}(undef, npairs)
    k = 1
    for e in 1:nelems, b in 1:4, a in 1:4
        I[k] = conn[a, e]
        J[k] = conn[b, e]
        k   += 1
    end

    # sparse() deduplicates and sorts; initialise nzval to 0 for assembly
    K = sparse(I, J, zeros(npairs), nnodes, nnodes)

    # --- Build elem_ptrs: map (a, b, e) → index into K.nzval ---
    #     CSC stores column b in nzval[colptr[b] : colptr[b+1]-1],
    #     with corresponding row indices in rowval[...].
    #     searchsortedfirst finds the position of row a in O(log nnz/col).
    elem_ptrs = Matrix{Int}(undef, 16, nelems)
    for e in 1:nelems
        for b in 1:4
            col   = conn[b, e]
            rng   = K.colptr[col] : K.colptr[col+1]-1
            for a in 1:4
                row = conn[a, e]
                pos = searchsortedfirst(K.rowval, row, first(rng), last(rng), Base.Forward)
                elem_ptrs[(b-1)*4 + a, e] = pos
            end
        end
    end

    K, elem_ptrs
end

# =============================================================================
# Step 2a — Scalar assembly
# =============================================================================
function assemble_scalar!(K, coords, conn, elem_ptrs, color_groups)
    fill!(K.nzval, 0.0)
    for group in color_groups
        for e in group
            Ke = element_Ke(coords, view(conn, :, e))
            @inbounds for b in 1:4, a in 1:4
                K.nzval[elem_ptrs[(b-1)*4+a, e]] += Ke[a, b]
            end
        end
    end
    K
end

# =============================================================================
# Step 2b — @threads assembly
#
#   Within a color group, no two elements share a node, so no two threads
#   write to the same nzval position.  The implicit barrier at the end of
#   @threads synchronises before the next color begins.
# =============================================================================
function assemble_threads!(K, coords, conn, elem_ptrs, color_groups)
    fill!(K.nzval, 0.0)
    for group in color_groups
        @threads for e in group
            Ke = element_Ke(coords, view(conn, :, e))
            @inbounds for b in 1:4, a in 1:4
                K.nzval[elem_ptrs[(b-1)*4+a, e]] += Ke[a, b]
            end
        end                             # ← implicit barrier
    end
    K
end

# =============================================================================
# Step 2c — @spawn / @sync assembly
#
#   Explicitly partition each color group into nthreads() chunks and spawn
#   one task per chunk.  @sync provides the barrier between colors.
#   Compared to @threads, this exposes work-stealing and lets tasks migrate
#   across threads; it is also composable with other async work.
# =============================================================================
function assemble_spawn!(K, coords, conn, elem_ptrs, color_groups)
    fill!(K.nzval, 0.0)
    nt = nthreads()
    for group in color_groups
        ne    = length(group)
        chunk = cld(ne, nt)
        @sync for t in 1:nt
            istart = (t-1)*chunk + 1
            iend   = min(t*chunk, ne)
            istart > ne && break
            @spawn begin
                @inbounds for idx in istart:iend
                    e = group[idx]
                    Ke = element_Ke(coords, view(conn, :, e))
                    for b in 1:4, a in 1:4
                        K.nzval[elem_ptrs[(b-1)*4+a, e]] += Ke[a, b]
                    end
                end
            end
        end                             # ← @sync barrier
    end
    K
end

# =============================================================================
# Driver
# =============================================================================
nx, ny = 512, 512

nnodes, coords, conn, colors = make_mesh(nx, ny)
color_groups = [findall(==(c), colors) for c in 1:4]

println("Building sparsity pattern and elem_ptrs …")
@time K, elem_ptrs = build_sparse_K(coords, conn, nnodes)
println("  K: $(nnodes)×$(nnodes), $(nnz(K)) non-zeros, $(size(conn,2)) elements")
println("  Threads: ", nthreads())
println()

# --- Verify all three results agree ---
K1 = copy(K);  assemble_scalar!(K1,  coords, conn, elem_ptrs, color_groups)
K2 = copy(K);  assemble_threads!(K2, coords, conn, elem_ptrs, color_groups)
K3 = copy(K);  assemble_spawn!(K3,   coords, conn, elem_ptrs, color_groups)

println("Max error (@threads) : ", maximum(abs.(K1.nzval .- K2.nzval)))
println("Max error (@spawn)   : ", maximum(abs.(K1.nzval .- K3.nzval)))

# Sanity check: row sums of assembled K should be zero on interior nodes
#   (Laplace with no BCs applied; each row of the stiffness matrix sums to 0)
row_sums = vec(sum(K1; dims=2))
interior  = [(j*(nx+1)+i+1) for j in 1:ny-1 for i in 1:nx-1]
println("Max |row-sum| on interior nodes: ", maximum(abs.(row_sums[interior])))
println()

# --- Benchmarks ---
print("Scalar   : "); @btime assemble_scalar!($K1, $coords, $conn, $elem_ptrs, $color_groups)
print("@threads : "); @btime assemble_threads!($K2, $coords, $conn, $elem_ptrs, $color_groups)
print("@spawn   : "); @btime assemble_spawn!($K3,  $coords, $conn, $elem_ptrs, $color_groups)

