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

# ============================================================================
# MFEM-Specific Functions
# ============================================================================

"""
Build CSR sparsity pattern for structured fine grid
ratio: refinement ratio (fine grid is (ratio+1) x (ratio+1) nodes)
"""
function build_fine_grid_sparsity(ratio::Int)
    numNodes = (ratio + 1) * (ratio + 1)

    matRowPtr = zeros(Int, numNodes + 1)
    matColIdx = Int[]

    for iy = 0:ratio
        for ix = 0:ratio
            nodeID = ix + iy * (ratio + 1) + 1  # 1-based indexing

            # Add neighbors in 9-point stencil
            if iy > 0
                if ix > 0
                    push!(matColIdx, nodeID - 1 - (ratio + 1))
                end
                push!(matColIdx, nodeID - (ratio + 1))
                if ix < ratio
                    push!(matColIdx, nodeID + 1 - (ratio + 1))
                end
            end
            if ix > 0
                push!(matColIdx, nodeID - 1)
            end
            push!(matColIdx, nodeID)
            if ix < ratio
                push!(matColIdx, nodeID + 1)
            end
            if iy < ratio
                if ix > 0
                    push!(matColIdx, nodeID - 1 + (ratio + 1))
                end
                push!(matColIdx, nodeID + (ratio + 1))
                if ix < ratio
                    push!(matColIdx, nodeID + 1 + (ratio + 1))
                end
            end

            matRowPtr[nodeID + 1] = length(matColIdx)
        end
    end

    # Convert to 1-based CSR
    matRowPtr[1] = 1
    for i = 2:numNodes+1
        matRowPtr[i] += 1
    end

    return matRowPtr, matColIdx
end

"""
Assemble fine grid matrix and RHS for a single coarse element
Uses optimized assemble_element! from FEMBase with workspace
"""
function assemble_fine_grid!(matValues::Vector{Float64}, rhs::Vector{Float64},
                            matRowPtr::Vector{Int}, matColIdx::Vector{Int},
                            elem::Q1Element{Dim, NNodes}, workspace::ElementWorkspace{Float64, Dim, NNodes},
                            ratio::Int,
                            x_corner::Vector{Float64}, y_corner::Vector{Float64},
                            ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    fill!(matValues, 0.0)
    fill!(rhs, 0.0)

    # Compute fine grid spacing
    hx = (x_corner[2] - x_corner[1]) / ratio
    hy = (y_corner[4] - y_corner[1]) / ratio

    # Assemble each fine element (reuse matrices across iterations)
    Ke_fine = zeros(4, 4)
    fe_fine = zeros(4)

    for iy = 0:ratio-1
        for ix = 0:ratio-1
            # Fine element corners
            x_fine = zeros(4)
            y_fine = zeros(4)
            x_fine[1] = x_corner[1] + ix * hx
            y_fine[1] = y_corner[1] + iy * hy
            x_fine[2] = x_fine[1] + hx
            y_fine[2] = y_fine[1]
            x_fine[3] = x_fine[2]
            y_fine[3] = y_fine[2] + hy
            x_fine[4] = x_fine[1]
            y_fine[4] = y_fine[3]

            # Assemble fine element using optimized kernel from FEMBase
            assemble_element!(Ke_fine, fe_fine, elem, workspace, x_fine, y_fine, ax_func, ay_func, f_func)

            # Fine element node list (1-based)
            nodeList = [
                ix + iy * (ratio + 1) + 1,
                ix + 1 + iy * (ratio + 1) + 1,
                ix + 1 + (iy + 1) * (ratio + 1) + 1,
                ix + (iy + 1) * (ratio + 1) + 1
            ]

            # Scatter to global fine grid arrays
            for i = 1:4
                gi = nodeList[i]
                rhs[gi] += fe_fine[i]

                for j = 1:4
                    gj = nodeList[j]
                    # Find position in sparse matrix
                    for k = matRowPtr[gi]:matRowPtr[gi + 1] - 1
                        if matColIdx[k] == gj
                            matValues[k] += Ke_fine[i, j]
                            break
                        end
                    end
                end
            end
        end
    end
end

"""
Compute MFEM basis functions and coarse element matrices
Returns: phi (basis functions on fine grid)
Uses optimized assembly with workspace
"""
function compute_mfem_element!(Ke_coarse::Matrix{Float64}, fe_coarse::Vector{Float64},
                              elem::Q1Element{Dim, NNodes}, workspace::ElementWorkspace{Float64, Dim, NNodes},
                              ratio::Int,
                              x_corner::Vector{Float64}, y_corner::Vector{Float64},
                              ax_func::Function, ay_func::Function, f_func::Function) where {Dim, NNodes}
    numNodes = (ratio + 1) * (ratio + 1)
    numVectors = 4  # Number of coarse basis functions
    numVectorsToSolve = 3  # Solve only 3, get 4th from partition of unity

    # Build fine grid sparsity pattern
    matRowPtr, matColIdx = build_fine_grid_sparsity(ratio)
    nnz = length(matColIdx)
    matValues = zeros(nnz)
    rhs_fine = zeros(numNodes)

    # Assemble fine grid matrix and RHS using optimized kernel
    assemble_fine_grid!(matValues, rhs_fine, matRowPtr, matColIdx, elem, workspace, ratio,
                       x_corner, y_corner, ax_func, ay_func, f_func)

    # Build DOF mapping (interior = free, boundary = fixed)
    globalToFree = fill(-1, numNodes)
    freeToGlobal = Int[]

    for iy = 0:ratio
        for ix = 0:ratio
            # Skip boundary nodes
            if ix == 0 || ix == ratio || iy == 0 || iy == ratio
                continue
            end
            nodeID = ix + iy * (ratio + 1) + 1
            push!(freeToGlobal, nodeID)
            globalToFree[nodeID] = length(freeToGlobal)
        end
    end

    nfree = length(freeToGlobal)

    # Initialize basis functions
    phi = zeros(numNodes, numVectors)

    # Set boundary conditions for first 3 basis functions (Q1 interpolation)
    for iy = 0:ratio
        eta = Float64(iy) / ratio

        # Left edge (x=0)
        ix = 0
        nodeID = ix + iy * (ratio + 1) + 1
        phi[nodeID, 1] = 1.0 - eta  # Basis 1: bottom-left corner
        phi[nodeID, 4] = eta        # Basis 4: top-left corner

        # Right edge (x=1)
        ix = ratio
        nodeID = ix + iy * (ratio + 1) + 1
        phi[nodeID, 2] = 1.0 - eta  # Basis 2: bottom-right corner
        phi[nodeID, 3] = eta        # Basis 3: top-right corner
    end

    for ix = 0:ratio
        xi = Float64(ix) / ratio

        # Bottom edge (y=0)
        iy = 0
        nodeID = ix + iy * (ratio + 1) + 1
        phi[nodeID, 1] = 1.0 - xi  # Basis 1: bottom-left corner
        phi[nodeID, 2] = xi        # Basis 2: bottom-right corner

        # Top edge (y=1)
        iy = ratio
        nodeID = ix + iy * (ratio + 1) + 1
        phi[nodeID, 3] = xi        # Basis 3: top-right corner
        phi[nodeID, 4] = 1.0 - xi  # Basis 4: top-left corner
    end

    # Build reduced system (free DOFs only)
    reduced_row_ptr = zeros(Int, nfree + 1)
    reduced_row_ptr[1] = 1

    for i = 1:nfree
        gi = freeToGlobal[i]
        count = 0
        for k = matRowPtr[gi]:matRowPtr[gi + 1] - 1
            gj = matColIdx[k]
            if globalToFree[gj] != -1
                count += 1
            end
        end
        reduced_row_ptr[i + 1] = reduced_row_ptr[i] + count
    end

    reduced_nnz = reduced_row_ptr[nfree + 1] - 1
    reduced_col_idx = zeros(Int, reduced_nnz)
    reduced_values = zeros(reduced_nnz)

    pos = 1
    for i = 1:nfree
        gi = freeToGlobal[i]
        for k = matRowPtr[gi]:matRowPtr[gi + 1] - 1
            gj = matColIdx[k]
            if globalToFree[gj] != -1
                reduced_col_idx[pos] = globalToFree[gj]
                reduced_values[pos] = matValues[k]
                pos += 1
            end
        end
    end

    # Build RHS for basis function solves
    btmp = zeros(nfree, numVectorsToSolve)

    for i = 1:nfree
        gi = freeToGlobal[i]
        for ir = 1:numVectorsToSolve
            sum_val = 0.0
            for k = matRowPtr[gi]:matRowPtr[gi + 1] - 1
                gj = matColIdx[k]
                sum_val += matValues[k] * phi[gj, ir]
            end
            btmp[i, ir] = -sum_val
        end
    end

    # Solve for interior basis functions
    K_reduced = SparseMatrixCSC(nfree, nfree, reduced_row_ptr, reduced_col_idx, reduced_values)
    utmp = K_reduced \ btmp

    # Copy solutions to phi
    for i = 1:nfree
        gi = freeToGlobal[i]
        for ir = 1:numVectorsToSolve
            phi[gi, ir] = utmp[i, ir]
        end
    end

    # Compute 4th basis function using partition of unity
    for i = 1:numNodes
        phi[i, 4] = 1.0 - phi[i, 1] - phi[i, 2] - phi[i, 3]
    end

    # Compute coarse element RHS: fe_coarse[ir] = phi[:,ir]' * rhs_fine
    for ir = 1:numVectors
        sum_val = 0.0
        for i = 1:numNodes
            sum_val += phi[i, ir] * rhs_fine[i]
        end
        fe_coarse[ir] = sum_val
    end

    # Compute coarse element stiffness: Ke_coarse[ir,jr] = phi[:,ir]' * K_fine * phi[:,jr]
    Kphi = zeros(numNodes, numVectors)

    # Matrix-vector products: Kphi = K_fine * phi
    for ir = 1:numVectors
        for i = 1:numNodes
            sum_val = 0.0
            for k = matRowPtr[i]:matRowPtr[i + 1] - 1
                j = matColIdx[k]
                sum_val += matValues[k] * phi[j, ir]
            end
            Kphi[i, ir] = sum_val
        end
    end

    # Compute Ke_coarse = phi' * Kphi
    for ir = 1:numVectors
        for jr = 1:numVectors
            sum_val = 0.0
            for i = 1:numNodes
                sum_val += phi[i, ir] * Kphi[i, jr]
            end
            Ke_coarse[ir, jr] = sum_val
        end
    end

    return phi  # Return basis functions for fine-scale reconstruction
end

# ============================================================================
# System Assembly with MFEM
# ============================================================================

"""
Assemble global system with MFEM elements (thread-parallel with coloring)
Returns: mat_values, rhs, basis_functions
Uses optimized kernels from FEMBase with thread-local workspaces
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

    # Thread-local element matrices and workspaces
    max_tid = Threads.maxthreadid()
    Ke_local = [zeros(4, 4) for _ = 1:max_tid]
    fe_local = [zeros(4) for _ = 1:max_tid]
    workspace_local = [ElementWorkspace{Float64, Dim, NNodes}() for _ = 1:max_tid]

    # Loop over colors
    for (ic, elements) in enumerate(e2e_colors)
        # All elements in this color can be assembled in parallel (no conflicts)
        Threads.@threads for iel in elements
            tid = Threads.threadid()
            @inbounds Ke = Ke_local[tid]
            @inbounds fe = fe_local[tid]
            @inbounds workspace = workspace_local[tid]

            # Get element nodes
            @inbounds nodes = mesh.cell_to_node[iel, :]
            @inbounds x = mesh.vertex_x[nodes]
            @inbounds y = mesh.vertex_y[nodes]

            # Assemble MFEM element and get basis functions using optimized kernel
            phi = compute_mfem_element!(Ke, fe, elem, workspace, ratio, x, y, ax_func, ay_func, f_func)
            basis_functions[iel] = phi

            # Scatter to global arrays (no race condition within same color)
            for i = 1:4
                @inbounds gi = nodes[i]
                @inbounds rhs[gi] += fe[i]
            end

            for i = 1:4
                @inbounds gi = nodes[i]
                for j = 1:4
                    @inbounds gj = nodes[j]
                    # Use binary search from FEMBase (assumes sorted columns)
                    k = find_matrix_position(n2n_row_ptr, n2n_col_idx, gi, gj)
                    @inbounds mat_values[k] += Ke[i, j]
                end
            end
        end
    end

    return mat_values, rhs, basis_functions
end

# ============================================================================
# Fine-Scale Solution Reconstruction
# ============================================================================

"""
Reconstruct fine-scale solution from coarse solution using MFEM basis functions
"""
function reconstruct_fine_solution(mesh::Mesh, coarse_solution::Vector{Float64},
                                  basis_functions::Vector{Matrix{Float64}}, ratio::Int)
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

    # Build fine mesh
    idx = 1
    for j = 0:ny_fine
        for i = 0:nx_fine
            fine_x[idx] = Float64(i) / nx_fine
            fine_y[idx] = Float64(j) / ny_fine
            idx += 1
        end
    end

    # Reconstruct solution element by element
    for iel = 1:nel
        # Coarse element indices
        iel_1d = iel - 1
        ix_coarse = iel_1d % nx_coarse
        iy_coarse = iel_1d ÷ nx_coarse

        # Coarse element nodes and solution values
        nodes = mesh.cell_to_node[iel, :]
        u_coarse_elem = coarse_solution[nodes]

        # Basis functions for this element
        phi = basis_functions[iel]

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
                for i = 1:4
                    u_val += phi[local_node, i] * u_coarse_elem[i]
                end
                fine_u[fine_node] = u_val
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

    # Thread info
    println("Number of threads: ", Threads.nthreads())
    println()

    # ========================================================================
    # Problem setup
    # ========================================================================

    # Material coefficients and forcing term (varying coefficients)
    ax(x, y, z) = 1.0
    ay(x, y, z) = 1.0
    function f(x, y, z)
        # Manufactured solution: u = sin(π*x) * sin(π*y)
        return 2.0 * π^2 * sin(π * x) * sin(π * y)
    end

    # ========================================================================
    # Mesh generation
    # ========================================================================

    nx, ny = 16, 16  # Coarse mesh (MFEM)
    ratio = 8        # Fine grid refinement ratio per coarse element

    println("Generating coarse mesh: $nx x $ny Q1 elements")
    println("MFEM ratio: $ratio (each coarse element has $(ratio)x$(ratio) fine elements)")
    println("Effective fine resolution: $(nx*ratio) x $(ny*ratio)")
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

    t0 = time()
    mat_values, rhs, basis_functions = assemble_system_mfem(mesh, elem, ratio, n2n_row_ptr, n2n_col_idx, e2e_colors,
                                                            ax, ay, f)
    assembly_time = time() - t0

    println("MFEM assembly complete")
    println("  Assembly time: ", @sprintf("%.2f ms", assembly_time * 1000))
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
        coarse_solution[free_to_global[i]] = sol_free[i]
    end

    t0 = time()
    fine_x, fine_y, fine_u = reconstruct_fine_solution(mesh, coarse_solution, basis_functions, ratio)
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
    println("  MFEM assembly:      ", @sprintf("%8.2f ms", assembly_time * 1000))
    println("  Boundary conditions:", @sprintf("%8.2f ms", bc_time * 1000))
    println("  Solve:              ", @sprintf("%8.2f ms", solve_time * 1000))
    println("  Reconstruction:     ", @sprintf("%8.2f ms", reconstruct_time * 1000))
    println("  " * "-"^68)
    total_time = mesh_time + conn_time + assembly_time + bc_time + solve_time + reconstruct_time
    println("  Total:              ", @sprintf("%8.2f ms", total_time * 1000))
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
