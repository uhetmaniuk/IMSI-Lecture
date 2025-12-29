# Examples

This directory contains standalone examples demonstrating various Kokkos programming patterns and FEM assembly techniques.

## Building

Examples are built as part of the main CMake build. From the repository root:

```bash
mkdir build && cd build
cmake .. -DKokkos_ENABLE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

Some examples require additional CMake options (see individual example descriptions).

## Examples

### Laplace1D.exe

**Source:** `Laplace1D.cpp`
**Requires:** `-DUSE_TACHO=ON`

Demonstrates solving a 1D Laplace equation using the Tacho sparse direct solver. Assembles a tridiagonal matrix and solves the linear system with timing benchmarks.

```bash
./examples/Laplace1D.exe
```

### simdPolar.exe

**Source:** `simdPolar.cpp`

Benchmarks different SIMD and parallelization strategies for a simple computation (`r = x*x + y*y + z*z`):
- Scalar loop (baseline)
- Compiler auto-vectorization (LLVM pragmas)
- Kokkos SIMD (`Kokkos::Experimental::native_simd`)
- Kokkos `parallel_for` with Serial and OpenMP execution spaces
- Kokkos `TeamPolicy` with `ThreadVectorRange`

Useful for understanding SIMD width on your hardware and comparing different parallelization approaches.

```bash
./examples/simdPolar.exe
```

### TeamPolicy.exe

**Source:** `team_policy_solution.cpp`

Demonstrates Kokkos `TeamPolicy` for hierarchical parallelism with both OpenMP and CUDA backends. Shows:
- League/team size configuration
- Nested `TeamThreadRange` parallelism
- Functor patterns for avoiding nested lambdas (required for CUDA)
- Atomic operations across threads

```bash
./examples/TeamPolicy.exe --kokkos-num-threads=8
```

### openmp_assembly.exe

**Source:** `openmp_assembly.cpp`
**Requires:** `-DKokkos_ENABLE_OPENMP=ON`

MsFEM/FEM assembly example using `ScaledLaplacian` with OpenMP execution space. Solves `-div(alpha * grad(u)) = f` on a 2D rectangular domain with homogeneous Dirichlet boundary conditions.

Demonstrates the portable assembly approach that works with any Kokkos execution space.

```bash
./examples/openmp_assembly.exe --kokkos-num-threads=8
```

### cuda_assembly.exe

**Source:** `cuda_assembly.cpp`
**Requires:** `-DKokkos_ENABLE_CUDA=ON`

CUDA-accelerated version of the MsFEM/FEM assembly example. Same problem as `openmp_assembly.exe` but runs on GPU.

```bash
./examples/cuda_assembly.exe
```

### Solver_Comparison.exe

**Source:** `Solver_Comparison.cpp`

Compares different sparse linear solver implementations and strategies.

```bash
./examples/Solver_Comparison.exe
```

## Configuration

The assembly examples (`openmp_assembly.cpp`, `cuda_assembly.cpp`) use `config_assembly.h` for shared configuration:
- `Scalar` type definition (float/double)
- Diffusion coefficient functors (`AxFunctor`, `AyFunctor`)
- Right-hand side functor (`FFunctor`)
