# Julia FEM Implementation

Julia implementations of the FEM assembly for comparison with the C++/Kokkos version.

## Files

### openmp_assembly.jl

Standard Q1 finite element assembly with thread-parallel execution.

Features:
- Bilinear Q1 shape functions with 2x2 Gauss quadrature
- Graph coloring for race-free parallel assembly
- Homogeneous Dirichlet boundary conditions
- Julia's built-in sparse direct solver

```bash
julia --threads=8 julia/openmp_assembly.jl
```

Output: `julia/solution.txt`

### mfem_assembly.jl

Multiscale Finite Element Method (MFEM) implementation.

Features:
- Coarse elements subdivided into fine grids (configurable ratio)
- Local basis functions via harmonic extension
- Partition of unity constraint
- Thread-parallel assembly with graph coloring

```bash
julia --threads=8 julia/mfem_assembly.jl
```

Output: `julia/mfem_solution.txt`

### mfem_cuda.jl

CUDA-accelerated MFEM implementation using CUDA.jl.

Features:
- One CUDA block per coarse element
- Block-cooperative PCG solver for interior basis functions
- Hybrid shared/global memory approach

```bash
julia julia/mfem_cuda.jl
```

## Configuration

Edit parameters in `main()`:
- `nx, ny`: Mesh resolution
- `ratio`: Fine grid refinement (for MFEM)
- Coefficient functions: `ax`, `ay`, `f`
