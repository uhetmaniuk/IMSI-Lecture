# IMSI-Lecture

Finite element method (FEM) teaching codebase using Kokkos for performance portability.

## Quick Start

```bash
mkdir build && cd build
cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON -DCMAKE_BUILD_TYPE=Release
cmake --build .
./main.exe --kokkos-num-threads=8
```

## CMake Options

- `-DKokkos_ENABLE_OPENMP=ON`: Enable OpenMP backend (recommended)
- `-DKokkos_ENABLE_CUDA=ON`: Enable CUDA backend
- `-DUSE_TACHO=ON`: Enable Tacho sparse solver
- `-DUSE_MKL=ON`: Use Intel MKL for BLAS/LAPACK

## Runtime Flags

- `--simd`: Enable SIMD vectorization for assembly
- `--nocolor`: Disable graph coloring (use atomics instead)
- `--kokkos-num-threads=N`: Set thread count
