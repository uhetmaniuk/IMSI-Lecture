# Source Code

## Core Components

| File | Description |
|------|-------------|
| `Mesh.h/cpp` | Mesh data structure (vertices, connectivity, boundary nodes) |
| `MeshUtils.h/cpp` | Mesh generation, graph operations (transpose, combine, coloring) |
| `Element.h` | Element type enumeration |
| `QuadratureRule.h/cpp` | Numerical integration rules |
| `Utils.h/cpp` | I/O utilities (GMSH output) |

## Finite Elements

| File | Description |
|------|-------------|
| `fe1DQ1.h`, `fe1DQ2.h` | 1D linear/quadratic elements |
| `fe2DQ1.h`, `fe2DQ2.h` | 2D bilinear/biquadratic quads |
| `fe3DQ1.h`, `fe3DQ2.h` | 3D trilinear/triquadratic hexahedra |

## Assembly

| File | Description |
|------|-------------|
| `ScaledLaplacian.h` | FEM assembly for `-div(alpha * grad(u)) = f` |
| `ScaledLaplacianHost.h` | Host-specific assembly routines |
| `FunctionExamples.h` | Example coefficient functions |

## Linear Algebra

| File | Description |
|------|-------------|
| `SparseMatrix.hpp` | CSR matrix operations |
| `SymmetricSparse.hpp` | Symmetric sparse LDL^T solver (wraps Fortran) |
| `PCG_Solver.h` | Preconditioned conjugate gradient solver |
| `SSOR_Solver.h` | SSOR preconditioner |
| `MathUtils.h` | Math utilities |
| `MathUtilsCuda.h` | CUDA-specific math utilities |

## Fortran Files (`*.f`)

Legacy sparse direct solver routines (SPOOLES-like). Called via `SymmetricSparse.hpp`.
