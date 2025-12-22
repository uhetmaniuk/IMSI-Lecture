#!/usr/bin/env julia
#
# @file mfem_assembly.jl
# @brief Julia implementation of 2D MFEM (Multiscale Finite Element Method)
#
# Solves: -∇·(α∇u) = f on a 2D rectangular domain
# with homogeneous Dirichlet boundary conditions
#
# Features:
# - 2D MFEM elements with local fine-grid resolution
# - Varying coefficients ax, ay, f
# - Thread-parallel assembly with graph coloring
# - Built-in sparse solver (\)
# - Static condensation for coarse-scale system
#

using SparseArrays
using LinearAlgebra
using Printf

# Fix for thread scaling: disable nested threading in BLAS
# This prevents oversubscription when using Threads.@threads
BLAS.set_num_threads(1)

# Load shared FEM base kernels
include("FEMBase.jl")
using .FEMBase

# Load PCG solver
include("PCG.jl")
using .PCG

# Load MFEM workspace (includes FEMBase and PCG)
include("MFEMWorkspace.jl")
using .MFEMWorkspace

# ============================================================================
# MFEM-Specific Functions
# ============================================================================

"""
Compute MFEM basis functions and coarse element matrices
Returns: phi (basis functions on fine grid), fine_asm_time, pcg_time, coarse_time
Uses optimized assembly with sparse K_ii (interior-interior) and K_b (boundary-all) matrices
Follows the C++ implementation in src/ScaledLaplacian.h (MFEMAssemblyFunctor::operator())

Note: Individual phase timings (fine_asm_time, pcg_time, coarse_time) are returned but not used
      in the current implementation. Instead, color-level wall-clock timing is used to match
      C++ methodology and avoid timing overcounting issues.
"""
function compute_mfem_element!(Ke_coarse::Matrix{Float64}, fe_coarse::Vector{Float64},
                              elem::Q1Element{Dim, NNodes}, workspace::MFEWorkspace{Float64, Dim, NNodes},
                              ratio::Int,
                              x_corner::Vector{Float64}, y_corner::Vector{Float64},
                              ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    numNodes = (ratio + 1) * (ratio + 1)
    numVectors = 4  # Number of coarse basis functions
    numVectorsToSolve = 3  # Solve only 3, get 4th from partition of unity

    # Build DOF mappings and sparsity patterns using reusable kernel from FEMBase
    nfree, nboundary, nnz_ii, nnz_b = build_mfem_local_sparsity!(workspace, ratio, ratio)

    # Extract arrays from workspace for later use
    globalToFree = workspace.globalToFree
    freeToGlobal = workspace.freeToGlobal
    globalToBoundary = workspace.globalToBoundary
    boundaryToGlobal = workspace.boundaryToGlobal
    colptr_ii = workspace.colptr_ii
    rowidx_ii = workspace.rowidx_ii
    nzval_ii = workspace.nzval_ii
    matRowPtr_b = workspace.matRowPtr_b
    matColIdx_b = workspace.matColIdx_b
    matValues_b = workspace.matValues_b

    # Use pre-allocated RHS and basis function arrays from workspace
    rhs_fine = workspace.rhs_fine
    btmp = workspace.btmp
    phi = workspace.phi
    utmp_init = workspace.utmp_init

    # Reset arrays
    fill!(rhs_fine, 0.0)
    fill!(btmp, 0.0)

    # Set corner values
    @inbounds phi[1, 1] = 1.0                                    # Bottom-left corner -> Basis 0
    @inbounds phi[ratio + 1, 2] = 1.0                           # Bottom-right corner -> Basis 1
    @inbounds phi[(ratio + 1) * (ratio + 1), 3] = 1.0           # Top-right corner -> Basis 2
    @inbounds phi[(ratio + 1) * ratio + 1, 4] = 1.0             # Top-left corner -> Basis 3

    # Set edge values (Q1 interpolation), skip corners
    for is = 1:ratio-1
        s = Float64(is) / ratio
        # Left edge (ix = 0)
        nodeID = is * (ratio + 1) + 1
        @inbounds phi[nodeID, 1] = 1.0 - s  # Basis 0
        @inbounds phi[nodeID, 4] = s        # Basis 3
        # Right edge (ix = ratio)
        nodeID = ratio + is * (ratio + 1) + 1
        @inbounds phi[nodeID, 2] = 1.0 - s  # Basis 1
        @inbounds phi[nodeID, 3] = s        # Basis 2
        # Bottom edge (iy = 0)
        nodeID = is + 1
        @inbounds phi[nodeID, 1] = 1.0 - s  # Basis 0
        @inbounds phi[nodeID, 2] = s        # Basis 1
        # Top edge (iy = ratio)
        nodeID = is + ratio * (ratio + 1) + 1
        @inbounds phi[nodeID, 4] = 1.0 - s  # Basis 3
        @inbounds phi[nodeID, 3] = s        # Basis 2
    end

    # ========================================================================
    # Phase 1: Fine-grid assembly
    # ========================================================================
    t_fine_start = time()

    # Assemble fine elements with scatter to K_ii, K_b, and btmp
    @inbounds hx = (x_corner[2] - x_corner[1]) / ratio
    @inbounds hy = (y_corner[4] - y_corner[1]) / ratio

    # Use pre-allocated arrays from element workspace (no redundant allocation!)
    x_fine = workspace.element.x_coords
    y_fine = workspace.element.y_coords
    nodeList = workspace.element.nodes

    for iy = 0:ratio-1
        for ix = 0:ratio-1
            # Fine element corners
            @inbounds x_fine[1] = x_corner[1] + ix * hx
            @inbounds y_fine[1] = y_corner[1] + iy * hy
            @inbounds x_fine[2] = x_fine[1] + hx
            @inbounds y_fine[2] = y_fine[1]
            @inbounds x_fine[3] = x_fine[2]
            @inbounds y_fine[3] = y_fine[2] + hy
            @inbounds x_fine[4] = x_fine[1]
            @inbounds y_fine[4] = y_fine[3]

            # Assemble fine element (Ke and fe are in workspace.element)
            assemble_element!(elem, workspace.element, x_fine, y_fine, ax_func, ay_func, f_func)

            # Extract Ke and fe from element workspace
            Ke_fine = workspace.element.Ke
            fe_fine = workspace.element.fe

            # Fine element node list
            nodeList[1] = ix + iy * (ratio + 1) + 1
            nodeList[2] = ix + 1 + iy * (ratio + 1) + 1
            nodeList[3] = ix + 1 + (iy + 1) * (ratio + 1) + 1
            nodeList[4] = ix + (iy + 1) * (ratio + 1) + 1

            # Scatter RHS
            for i = 1:4
                @inbounds gi = nodeList[i]
                @inbounds rhs_fine[gi] += fe_fine[i]
            end

            # Scatter stiffness to K_ii, K_b, and btmp (hoisted branches for performance)
            for i = 1:4
                @inbounds iGlobal = nodeList[i]
                @inbounds iFree = globalToFree[iGlobal]

                if iFree != -1
                    # Interior row - process all j together
                    for j = 1:4
                        @inbounds jGlobal = nodeList[j]
                        @inbounds jFree = globalToFree[jGlobal]
                        @inbounds k_val = Ke_fine[i, j]

                        if jFree != -1
                            # Interior-interior: add to K_ii
                            k = find_matrix_position(colptr_ii, rowidx_ii, iFree, jFree)
                            @inbounds nzval_ii[k] += k_val
                        else
                            # Interior-boundary: add to btmp (RHS for static condensation)
                            @simd for ir = 1:numVectorsToSolve
                                @inbounds phi_val = phi[jGlobal, ir]
                                @inbounds btmp[iFree, ir] -= k_val * phi_val
                            end
                        end
                    end
                else
                    # Boundary row - process all j together
                    @inbounds iBoundary = globalToBoundary[iGlobal]
                    for j = 1:4
                        @inbounds jGlobal = nodeList[j]
                        @inbounds k_val = Ke_fine[i, j]
                        k = find_matrix_position(matRowPtr_b, matColIdx_b, iBoundary, jGlobal)
                        @inbounds matValues_b[k] += k_val
                    end
                end
            end
        end
    end

    fine_asm_time = time() - t_fine_start

    # ========================================================================
    # Phase 2: Build initial guess and solve for interior basis functions
    # ========================================================================
    t_pcg_start = time()

    # Build initial guess for PCG based on Q1 shape functions
    for i = 1:nfree
        @inbounds gi = freeToGlobal[i]

        # Compute normalized coordinates (xi, eta) of this fine node within coarse element
        # Fine nodes are arranged in a (ratio+1) x (ratio+1) grid
        iy_fine = (gi - 1) ÷ (ratio + 1)
        ix_fine = (gi - 1) % (ratio + 1)
        xi = Float64(ix_fine) / ratio   # Normalized coordinate [0, 1]
        eta = Float64(iy_fine) / ratio  # Normalized coordinate [0, 1]

        # Evaluate Q1 shape functions at this position
        # Q1 nodes: [1] bottom-left, [2] bottom-right, [3] top-right, [4] top-left
        N1 = (1.0 - xi) * (1.0 - eta)  # Basis 1 (bottom-left)
        N2 = xi * (1.0 - eta)          # Basis 2 (bottom-right)
        N3 = xi * eta                  # Basis 3 (top-right)

        # Initial guess for interior values of first 3 basis functions
        @inbounds utmp_init[i, 1] = N1
        @inbounds utmp_init[i, 2] = N2
        @inbounds utmp_init[i, 3] = N3
    end

    # Solve for interior basis functions using in-place PCG with SSOR preconditioning
    # Use only the filled portions of the pre-allocated arrays
    K_ii = SparseMatrixCSC(nfree, nfree, colptr_ii[1:nfree+1], rowidx_ii[1:nnz_ii], nzval_ii[1:nnz_ii])

    # Copy initial guess to utmp (will be modified in-place)
    utmp = workspace.utmp
    @. utmp = utmp_init

    # Solve in-place (ZERO allocations!)
    pcg_info = pcg_solve_ssor!(utmp, workspace.pcg_workspace, K_ii, btmp,
                               omega=1.0, num_ssor_sweeps=1,
                               tol=1e-12, maxiter=1000, verbose=false)

    # Check convergence
    #if !pcg_info.converged
    #    @warn "PCG-SSOR did not converge: iterations=$(pcg_info.iterations), residual=$(pcg_info.residual_norm)"
    #end

    pcg_time = time() - t_pcg_start

    # ========================================================================
    # Phase 3: Compute coarse element matrices and RHS
    # ========================================================================
    t_coarse_start = time()

    # Reconstruct full basis functions: copy interior solutions and compute 4th basis
    for i = 1:nfree
        @inbounds gi = freeToGlobal[i]
        sum_val = 0.0
        for ir = 1:numVectorsToSolve
            @inbounds phi[gi, ir] = utmp[i, ir]
            sum_val += utmp[i, ir]
        end
        # 4th basis function from partition of unity
        @inbounds phi[gi, 4] = 1.0 - sum_val
    end

    # Compute coarse element RHS: fe_coarse[ir] = phi[:,ir]' * rhs_fine
    for ir = 1:numVectors
        sum_val = 0.0
        @simd for i = 1:numNodes
            @inbounds sum_val += phi[i, ir] * rhs_fine[i]
        end
        @inbounds fe_coarse[ir] = sum_val
    end

    # Compute coarse element stiffness using Schur complement: Ke_coarse = phi^T * K * phi
    # This is computed efficiently using K_b (boundary rows only):
    # kele = phi_b^T * K_b * phi
    # where K_b includes boundary-interior and boundary-boundary coupling
    fill!(Ke_coarse, 0.0)

    for iBoundary = 1:nboundary
        @inbounds iGlobal = boundaryToGlobal[iBoundary]
        @inbounds rowBegin = matRowPtr_b[iBoundary]
        @inbounds rowEnd = matRowPtr_b[iBoundary + 1] - 1

        for k = rowBegin:rowEnd
            @inbounds jGlobal = matColIdx_b[k]  # K_b uses global node numbering for columns
            @inbounds k_val = matValues_b[k]

            # Accumulate contribution to Ke_coarse for all basis function pairs
            # Hoist phi_i * k_val out of jr loop for better performance
            for ir = 1:numVectors
                @inbounds phi_i_k = phi[iGlobal, ir] * k_val
                @simd for jr = 1:numVectors
                    @inbounds phi_j = phi[jGlobal, jr]
                    @inbounds Ke_coarse[ir, jr] += phi_i_k * phi_j
                end
            end
        end
    end

    coarse_time = time() - t_coarse_start

    return phi, fine_asm_time, pcg_time, coarse_time
end

# ============================================================================
# System Assembly with MFEM
# ============================================================================

"""
Assemble global system with MFEM elements (thread-parallel with coloring)
Returns: mat_values, rhs, basis_functions, kernel_time, fine_time, pcg_time, coarse_time, scatter_time
Uses optimized kernels from FEMBase with thread-local workspaces

Timing methodology: For each color, finds the slowest thread and uses its phase breakdown.
This ensures phases sum correctly to the total time and provides detailed performance insights.
"""
function assemble_system_mfem(mesh::Mesh, elem::Q1Element{Dim, NNodes}, ratio::Int,
                             n2n_row_ptr::Vector{Int}, n2n_col_idx::Vector{Int},
                             e2e_colors::Vector{Vector{Int}},
                             ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    nnodes = length(mesh.vertex_x)
    nnz = length(n2n_col_idx)
    nel = size(mesh.cell_to_node, 1)

    # Initialize global arrays
    mat_values = zeros(nnz)
    rhs = zeros(nnodes)

    # Store basis functions for each element (for fine-scale reconstruction)
    numFineNodes = (ratio + 1) * (ratio + 1)
    basis_functions = [zeros(numFineNodes, 4) for _ = 1:nel]

    # Thread-local workspaces (Ke_coarse and fe_coarse are now in workspace)
    max_tid = Threads.maxthreadid()
    workspace_local = [MFEWorkspace{Float64, Dim, NNodes}(ratio) for _ = 1:max_tid]

    # Timing accumulators (sum across colors)
    # Julia can provide detailed breakdown by tracking per-thread phase timings
    total_kernel_time = 0.0
    total_fine_time = 0.0
    total_pcg_time = 0.0
    total_coarse_time = 0.0
    total_scatter_time = 0.0

    # Loop over colors
    for (ic, elements) in enumerate(e2e_colors)
        # Time this color's assembly (wall-clock time)
        color_start = time()

        # Thread-local timing arrays for phase breakdown
        thread_total_times = zeros(max_tid)
        thread_fine_times = zeros(max_tid)
        thread_pcg_times = zeros(max_tid)
        thread_coarse_times = zeros(max_tid)

        # All elements in this color can be assembled in parallel (no conflicts)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            @inbounds workspace = workspace_local[tid]

            # Use pre-allocated coarse element arrays from workspace (avoid allocation!)
            nodes = workspace.nodes_coarse
            x = workspace.x_coarse
            y = workspace.y_coarse
            Ke = workspace.Ke_coarse
            fe = workspace.fe_coarse

            # Extract element data (manual copy, no allocation)
            @inbounds for i = 1:4
                nodes[i] = mesh.cell_to_node[iel, i]
            end
            @inbounds for i = 1:4
                ni = nodes[i]
                x[i] = mesh.vertex_x[ni]
                y[i] = mesh.vertex_y[ni]
            end

            # Assemble MFEM element and get basis functions using optimized kernel
            t_elem_start = time()
            @inbounds basis_functions[iel], t_fine, t_pcg, t_coarse = compute_mfem_element!(Ke, fe, elem, workspace, ratio, x, y, ax_func, ay_func, f_func)

            # Scatter to global arrays (no race condition within same color)
            for i = 1:4
                @inbounds gi = nodes[i]
                @inbounds rhs[gi] += fe[i]
                for j = 1:4
                    @inbounds gj = nodes[j]
                    # Use binary search from FEMBase (assumes sorted columns)
                    k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                    @inbounds mat_values[k] += Ke[i, j]
                end
            end
            t_elem_total = time() - t_elem_start

            # Accumulate timing for this thread
            thread_total_times[tid] += t_elem_total
            thread_fine_times[tid] += t_fine
            thread_pcg_times[tid] += t_pcg
            thread_coarse_times[tid] += t_coarse
        end

        # Wall-clock time for this color
        color_time = time() - color_start

        # Find the slowest thread (determines wall-clock time for this color)
        slowest_tid = argmax(thread_total_times)

        # Use the slowest thread's phase breakdown (ensures phases sum to total)
        color_fine = thread_fine_times[slowest_tid]
        color_pcg = thread_pcg_times[slowest_tid]
        color_coarse = thread_coarse_times[slowest_tid]
        color_compute = color_fine + color_pcg + color_coarse
        color_scatter = thread_total_times[slowest_tid] - color_compute

        # Sum across colors (colors are processed sequentially)
        total_kernel_time += color_time
        total_fine_time += color_fine
        total_pcg_time += color_pcg
        total_coarse_time += color_coarse
        total_scatter_time += color_scatter
    end

    return mat_values, rhs, basis_functions, total_kernel_time, total_fine_time, total_pcg_time, total_coarse_time, total_scatter_time
end

# ============================================================================
# Fine-Scale Solution Reconstruction
# ============================================================================

"""
Reconstruct fine-scale solution from coarse solution using MFEM basis functions
Uses graph coloring to avoid race conditions on shared fine nodes
"""
function reconstruct_fine_solution(mesh::Mesh, coarse_solution::Vector{Float64},
                                  basis_functions::Vector{Matrix{Float64}},
                                  e2e_colors::Vector{Vector{Int}}, ratio::Int)
    nel = size(mesh.cell_to_node, 1)
    nx_coarse = Int(sqrt(nel))  # Assuming square mesh
    ny_coarse = nx_coarse

    # Fine mesh dimensions
    nx_fine = nx_coarse * ratio
    ny_fine = ny_coarse * ratio
    nfine = (nx_fine + 1) * (ny_fine + 1)

    # Fine mesh coordinates and solution
    fine_x = zeros(nfine)
    fine_y = zeros(nfine)
    fine_u = zeros(nfine)

    # Build fine mesh (parallel - independent writes)
    Threads.@threads for j = 0:ny_fine
        for i = 0:nx_fine
            idx = (i + 1) + j * (nx_fine + 1)
            @inbounds fine_x[idx] = Float64(i) / nx_fine
            @inbounds fine_y[idx] = Float64(j) / ny_fine
        end
    end

    # Thread-local temporary arrays (avoid allocations in hot loop)
    max_tid = Threads.maxthreadid()
    nodes_local = [zeros(Int, 4) for _ = 1:max_tid]
    u_local = [zeros(4) for _ = 1:max_tid]

    # Reconstruct solution element by element (parallel with coloring to avoid races)
    for (ic, elements) in enumerate(e2e_colors)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            nodes = nodes_local[tid]
            u_coarse_elem = u_local[tid]

            # Coarse element indices
            iel_1d = iel - 1
            ix_coarse = iel_1d % nx_coarse
            iy_coarse = iel_1d ÷ nx_coarse

            # Extract coarse element nodes and solution values (no allocation!)
            @inbounds for i = 1:4
                nodes[i] = mesh.cell_to_node[iel, i]
                u_coarse_elem[i] = coarse_solution[nodes[i]]
            end

            # Basis functions for this element
            @inbounds phi = basis_functions[iel]

            # Reconstruct fine solution: u_fine = phi * u_coarse
            for iy_local = 0:ratio
                for ix_local = 0:ratio
                    # Fine node global indices
                    ix_fine = ix_coarse * ratio + ix_local
                    iy_fine = iy_coarse * ratio + iy_local
                    fine_node = ix_fine + iy_fine * (nx_fine + 1) + 1

                    # Local fine node index
                    local_node = ix_local + iy_local * (ratio + 1) + 1

                    # Interpolate: u_fine = sum_i phi[local_node, i] * u_coarse[i]
                    u_val = 0.0
                    @simd for i = 1:4
                        @inbounds u_val += phi[local_node, i] * u_coarse_elem[i]
                    end
                    @inbounds fine_u[fine_node] = u_val
                end
            end
        end
    end

    return fine_x, fine_y, fine_u
end

# ============================================================================
# Main Program
# ============================================================================

function main()
    println("="^70)
    println("Julia 2D MFEM Assembly (Thread-Parallel with Coloring)")
    println("="^70)
    println()
    println("Usage: julia -t <threads> mfem_assembly.jl [nx] [ny] [ratio]")
    println("  nx:     Coarse mesh elements in x-direction (default: 48)")
    println("  ny:     Coarse mesh elements in y-direction (default: nx)")
    println("  ratio:  Fine grid refinement ratio (default: 64)")
    println()

    # Thread info
    println("Number of threads: ", Threads.nthreads())
    println()

    # ========================================================================
    # Problem setup
    # ========================================================================

    # Material coefficients and forcing term
    # Set USE_VARYING_COEFFICIENTS to true to test with non-constant diffusion
    USE_VARYING_COEFFICIENTS = true

    # Define coefficient functions (use ternary operator to switch behavior)
    ax(x, y, z) = USE_VARYING_COEFFICIENTS ? (1.0 + 100 * cos(150 * x) * cos(150 * x) * sin(150 * y) * sin(150 * y)) : 1.0
    ay(x, y, z) = USE_VARYING_COEFFICIENTS ? (1.0 + 100 * cos(150 * x) * cos(150 * x) * sin(150 * y) * sin(150 * y)) : 1.0

    # Right-hand side (same for both cases)
    function f(x, y, z)
        # Manufactured solution: u = sin(π*x) * sin(π*y)
        # Note: For varying coefficients, this RHS is no longer exact
        return USE_VARYING_COEFFICIENTS ? sin(x) * sin(y) : 2.0 * π^2 * sin(π * x) * sin(π * y)
    end

    # ========================================================================
    # Mesh generation
    # ========================================================================

    # Parse command-line arguments: julia mfem_assembly.jl [nx] [ny] [ratio]
    # Defaults: nx=48, ny=48, ratio=64
    nx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 48
    ny = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nx  # Default ny=nx if not specified
    ratio = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 64

    println("Generating coarse mesh: $nx x $ny MFEM elements")
    println("MFEM ratio: $ratio (each coarse element has $(ratio)x$(ratio) fine elements)")
    println("Effective fine resolution: $(nx*ratio) x $(ny*ratio)")
    if USE_VARYING_COEFFICIENTS
        println("Diffusion coefficients: VARYING (1 + 100 * cos(150*x)^2 * sin(150*y)^2 terms)")
    else
        println("Diffusion coefficients: CONSTANT (ax=ay=1)")
    end
    println()

    t0 = time()
    mesh = generate_mesh(nx, ny)
    mesh_time = time() - t0

    println("  Number of coarse elements: ", size(mesh.cell_to_node, 1))
    println("  Number of coarse nodes:    ", length(mesh.vertex_x))
    println("  Mesh generation:           ", @sprintf("%.2f ms", mesh_time * 1000))
    println()

    # ========================================================================
    # Build connectivity and coloring
    # ========================================================================

    # Warm-up run to eliminate JIT compilation overhead
    println("Running connectivity JIT warm-up...")
    t_warmup = time()
    _, _, _ = build_mesh_connectivity(mesh)
    warmup_time = time() - t_warmup
    println("  Warm-up time: ", @sprintf("%.2f ms", warmup_time * 1000))
    println()

    println("Building mesh connectivity...")
    t0 = time()
    n2n_row_ptr, n2n_col_idx, e2e_colors = build_mesh_connectivity(mesh)
    conn_time = time() - t0

    nnz = length(n2n_col_idx)
    println("  Matrix size: ", length(mesh.vertex_x), " x ", length(mesh.vertex_x))
    println("  Non-zeros:   ", nnz)
    println("  Number of colors: ", length(e2e_colors))
    for ic = 1:length(e2e_colors)
        println("    Color $ic: ", length(e2e_colors[ic]), " elements")
    end
    println("  Connectivity time: ", @sprintf("%.2f ms", conn_time * 1000))
    println()

    # ========================================================================
    # Assembly
    # ========================================================================

    println("="^70)
    println("Starting MFEM Assembly")
    println("="^70)

    elem = Q1Element()

    # Warm-up run to eliminate JIT compilation overhead
    println("Running JIT warm-up...")
    t_warmup = time()
    _, _, _, _, _, _, _, _ = assemble_system_mfem(mesh, elem, ratio, n2n_row_ptr, n2n_col_idx, e2e_colors,
                                                   ax, ay, f)
    warmup_time = time() - t_warmup
    println("  Warm-up time: ", @sprintf("%.2f ms", warmup_time * 1000))
    println()

    # Measure total assembly time (including all overhead)
    t0 = time()
    mat_values, rhs, basis_functions, kernel_time, fine_time, pcg_time, coarse_time, scatter_time =
        assemble_system_mfem(mesh, elem, ratio, n2n_row_ptr, n2n_col_idx, e2e_colors, ax, ay, f)
    total_time = time() - t0

    # Overhead = total_time - kernel_time (includes startup, memory allocation, etc.)
    overhead_time = total_time - kernel_time

    println("MFEM assembly complete")
    println("  Total assembly time:   ", @sprintf("%.2f ms", total_time * 1000))
    println("  ├─ Fine assembly:      ", @sprintf("%.2f ms (%.1f%%)", fine_time * 1000, 100*fine_time/total_time))
    println("  ├─ PCG solve:          ", @sprintf("%.2f ms (%.1f%%)", pcg_time * 1000, 100*pcg_time/total_time))
    println("  ├─ Coarse computation: ", @sprintf("%.2f ms (%.1f%%)", coarse_time * 1000, 100*coarse_time/total_time))
    println("  ├─ Global scatter:     ", @sprintf("%.2f ms (%.1f%%)", scatter_time * 1000, 100*scatter_time/total_time))
    println("  ├─ Kernel time:        ", @sprintf("%.2f ms (%.1f%%)", kernel_time * 1000, 100*kernel_time/total_time))
    println("  └─ Overhead:           ", @sprintf("%.2f ms (%.1f%%)", overhead_time * 1000, 100*overhead_time/total_time))
    println()

    # ========================================================================
    # Apply boundary conditions
    # ========================================================================

    println("Applying boundary conditions...")
    t0 = time()
    reduced_row_ptr, reduced_col_idx, reduced_values, reduced_rhs, free_to_global =
        apply_boundary_conditions(n2n_row_ptr, n2n_col_idx, mat_values, rhs, mesh.boundary_nodes)
    bc_time = time() - t0

    nfree = length(free_to_global)
    println("  Total DOFs:    ", length(mesh.vertex_x))
    println("  Boundary DOFs: ", length(mesh.boundary_nodes))
    println("  Free DOFs:     ", nfree)
    println("  Reduced nnz:   ", length(reduced_col_idx))
    println("  BC time:       ", @sprintf("%.2f ms", bc_time * 1000))
    println()

    # ========================================================================
    # Solve
    # ========================================================================

    println("Solving linear system...")

    # Build sparse matrix
    K = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)

    t0 = time()
    sol_free = K \ reduced_rhs
    solve_time = time() - t0

    println("  Solve time: ", @sprintf("%.2f ms", solve_time * 1000))
    println()

    # ========================================================================
    # Reconstruct Fine-Scale Solution
    # ========================================================================

    println("Reconstructing fine-scale solution...")

    # Expand coarse solution to full DOFs (BCs = 0)
    coarse_solution = zeros(length(mesh.vertex_x))
    for i = 1:nfree
        @inbounds coarse_solution[free_to_global[i]] = sol_free[i]
    end

    t0 = time()
    fine_x, fine_y, fine_u = reconstruct_fine_solution(mesh, coarse_solution, basis_functions, e2e_colors, ratio)
    reconstruct_time = time() - t0

    println("  Fine mesh size: ", length(fine_u), " nodes")
    println("  Reconstruction time: ", @sprintf("%.2f ms", reconstruct_time * 1000))
    println()

    # ========================================================================
    # Output
    # ========================================================================

    println("Writing solutions to files...")

    # Write coarse solution
    write_solution("mfem_solution_coarse.txt", mesh, coarse_solution)
    println("  Coarse solution written to: mfem_solution_coarse.txt")

    # Write fine solution
    open("mfem_solution_fine.txt", "w") do io
        println(io, "# x y u")
        for i = 1:length(fine_u)
            @printf(io, "%.6f %.6f %.6e\n", fine_x[i], fine_y[i], fine_u[i])
        end
    end
    println("  Fine solution written to: mfem_solution_fine.txt")
    println()

    # ========================================================================
    # Performance summary
    # ========================================================================

    println("="^70)
    println("Performance Summary")
    println("="^70)
    println("  Mesh generation:    ", @sprintf("%8.2f ms", mesh_time * 1000))
    println("  Connectivity/color: ", @sprintf("%8.2f ms", conn_time * 1000))
    println("  MFEM assembly:      ", @sprintf("%8.2f ms", total_time * 1000))
    println("    ├─ Fine assembly: ", @sprintf("%8.2f ms (%.1f%%)", fine_time * 1000, 100*fine_time/total_time))
    println("    ├─ PCG solve:     ", @sprintf("%8.2f ms (%.1f%%)", pcg_time * 1000, 100*pcg_time/total_time))
    println("    ├─ Coarse compute:", @sprintf("%8.2f ms (%.1f%%)", coarse_time * 1000, 100*coarse_time/total_time))
    println("    ├─ Scatter:       ", @sprintf("%8.2f ms (%.1f%%)", scatter_time * 1000, 100*scatter_time/total_time))
    println("    ├─ Kernel total:  ", @sprintf("%8.2f ms (%.1f%%)", kernel_time * 1000, 100*kernel_time/total_time))
    println("    └─ Overhead:      ", @sprintf("%8.2f ms (%.1f%%)", overhead_time * 1000, 100*overhead_time/total_time))
    println("  Boundary conditions:", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Solve:              ", @sprintf("%8.2f ms", solve_time * 1000))
    println("  Reconstruction:     ", @sprintf("%8.2f ms", reconstruct_time * 1000))
    println("  " * "-"^68)
    total_workflow_time = mesh_time + conn_time + total_time + bc_time + solve_time + reconstruct_time
    println("  Total:              ", @sprintf("%8.2f ms", total_workflow_time * 1000))
    println()

    # Additional metrics for scaling studies
    println("="^70)
    println("Metrics for Scaling Studies")
    println("="^70)
    println("  Problem size:")
    println("    Coarse mesh:      $nx × $ny = ", nx*ny, " elements")
    println("    MFEM ratio:       $ratio")
    println("    Effective fine:   $(nx*ratio) × $(ny*ratio) = ", (nx*ratio)*(ny*ratio), " elements")
    println("    DOFs (coarse):    ", length(mesh.vertex_x))
    println("    Free DOFs:        ", nfree)
    println()
    println("  Assembly breakdown (detailed):")
    println("    Fine assembly:    ", @sprintf("%8.2f ms (%.1f%%)", fine_time * 1000, 100*fine_time/total_time))
    println("    PCG solve:        ", @sprintf("%8.2f ms (%.1f%%)", pcg_time * 1000, 100*pcg_time/total_time))
    println("    Coarse compute:   ", @sprintf("%8.2f ms (%.1f%%)", coarse_time * 1000, 100*coarse_time/total_time))
    println("    Scatter:          ", @sprintf("%8.2f ms (%.1f%%)", scatter_time * 1000, 100*scatter_time/total_time))
    println("    Kernel total:     ", @sprintf("%8.2f ms (%.1f%%)", kernel_time * 1000, 100*kernel_time/total_time))
    println("    Overhead:         ", @sprintf("%8.2f ms (%.1f%%)", overhead_time * 1000, 100*overhead_time/total_time))
    println("    Total assembly:   ", @sprintf("%8.2f ms", total_time * 1000))
    println()
    println("  Parallel efficiency:")
    println("    Threads used:     ", Threads.nthreads())
    println("    Number of colors: ", length(e2e_colors))
    nel_total = size(mesh.cell_to_node, 1)
    println("    Elements/thread:  ", @sprintf("%.1f", nel_total / Threads.nthreads()))
    println()

    # Solution statistics
    println("Coarse solution statistics:")
    minval = minimum(sol_free)
    maxval = maximum(sol_free)
    avgval = sum(sol_free) / length(sol_free)
    println("  Min value: ", @sprintf("%.6e", minval))
    println("  Max value: ", @sprintf("%.6e", maxval))
    println("  Avg value: ", @sprintf("%.6e", avgval))
    println()

    println("Fine solution statistics:")
    minval_fine = minimum(fine_u)
    maxval_fine = maximum(fine_u)
    avgval_fine = sum(fine_u) / length(fine_u)
    println("  Min value: ", @sprintf("%.6e", minval_fine))
    println("  Max value: ", @sprintf("%.6e", maxval_fine))
    println("  Avg value: ", @sprintf("%.6e", avgval_fine))
    println()
end

# Run main program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
