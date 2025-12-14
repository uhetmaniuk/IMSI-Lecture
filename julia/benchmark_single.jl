include("openmp_assembly.jl")
using Printf

# Problem setup
ax(x, y, z) = 1.0 + 0.5 * sin(2π * x) * cos(2π * y)
ay(x, y, z) = 1.0 + 0.5 * cos(2π * x) * sin(2π * y)
f(x, y, z) = 2.0 * π^2 * sin(π * x) * sin(π * y)

# Mesh generation
nx, ny = 512, 512
mesh = generate_mesh(nx, ny)

# Warmup
n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
elem = Q1Element()
mat_values, rhs, color_times, lookup_time = assemble_system(
    mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f
)

# Actual timing runs (3 trials)
println("Threads: ", Threads.nthreads())
mesh_times = Float64[]
conn_times = Float64[]
lookup_times = Float64[]
assembly_times = Float64[]
bc_times = Float64[]
solve_times = Float64[]

for trial = 1:3
    # Connectivity
    t0 = time()
    n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
    push!(conn_times, (time() - t0) * 1000)

    # Assembly
    t0 = time()
    mat_values, rhs, color_times, lookup_time = assemble_system(
        mesh, elem, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f
    )
    assembly_time = (time() - t0) * 1000
    push!(assembly_times, assembly_time)
    push!(lookup_times, lookup_time * 1000)

    # Boundary conditions
    t0 = time()
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global =
        apply_boundary_conditions(n2n_row_ptr, n2n_col_idx, mat_values, rhs, mesh.boundary_nodes)
    push!(bc_times, (time() - t0) * 1000)

    # Solve
    nfree = length(free_to_global)
    K = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)
    t0 = time()
    sol_free = K \ reduced_rhs
    push!(solve_times, (time() - t0) * 1000)
end

# Report median times
using Statistics
println(@sprintf("Connectivity: %.2f ms (median of 3)", median(conn_times)))
println(@sprintf("Lookup table: %.2f ms (median of 3)", median(lookup_times)))
println(@sprintf("Assembly:     %.2f ms (median of 3)", median(assembly_times)))
println(@sprintf("BC:           %.2f ms (median of 3)", median(bc_times)))
println(@sprintf("Solve:        %.2f ms (median of 3)", median(solve_times)))
total = median(conn_times) + median(assembly_times) + median(bc_times) + median(solve_times)
println(@sprintf("Total:        %.2f ms", total))
