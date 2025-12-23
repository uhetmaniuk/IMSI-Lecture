# C++ CUDA Scaling Study Guide

This guide explains how to run comprehensive CUDA scaling studies for different MFEM ratios.

## Overview

Four separate scripts are provided for testing different MFEM ratios:
- `run_cpp_cuda_r8.sh` - Ratio 8 (fine mesh: 8×8 per coarse element)
- `run_cpp_cuda_r16.sh` - Ratio 16 (fine mesh: 16×16 per coarse element)
- `run_cpp_cuda_r32.sh` - Ratio 32 (fine mesh: 32×32 per coarse element)
- `run_cpp_cuda_r64.sh` - Ratio 64 (fine mesh: 64×64 per coarse element)

Each script tests mesh sizes from 16×16 to 512×512 and compares CUDA vs OpenMP performance.

## IMPORTANT: Recompilation Required

The MFEM ratio is **hardcoded** in the source code. You must recompile for each ratio.

### Steps for Each Ratio

#### 1. Edit Source Code

Edit `src/ScaledLaplacian.h` line 144:

For ratio=8:
```cpp
static constexpr int ratio = 8;
```

For ratio=16:
```cpp
static constexpr int ratio = 16;
```

For ratio=32:
```cpp
static constexpr int ratio = 32;
```

For ratio=64:
```cpp
static constexpr int ratio = 64;
```

#### 2. Recompile

```bash
cd build
cmake .. -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON  # or your GPU architecture
make cuda_assembly openmp_assembly -j8
cd ..
```

**GPU Architecture Options:**
- `-DKokkos_ARCH_AMPERE80=ON` for NVIDIA A100
- `-DKokkos_ARCH_HOPPER90=ON` for NVIDIA H100
- `-DKokkos_ARCH_VOLTA70=ON` for NVIDIA V100
- See Kokkos documentation for other GPUs

#### 3. Run Script

```bash
./run_cpp_cuda_r8.sh   # For ratio=8
```

The script will:
1. Prompt you to confirm recompilation
2. Run CUDA tests (16×16 to 512×512 meshes)
3. Run OpenMP baseline tests
4. Generate summary report

#### 4. Repeat for Other Ratios

For each ratio (8, 16, 32, 64):
- Edit `src/ScaledLaplacian.h` (change ratio value)
- Recompile (`cd build && make cuda_assembly openmp_assembly`)
- Run corresponding script

## Output Directories

Results are saved to separate directories:
- `cpp_cuda_scaling_results_r8/` - Ratio 8 results
- `cpp_cuda_scaling_results_r16/` - Ratio 16 results
- `cpp_cuda_scaling_results_r32/` - Ratio 32 results
- `cpp_cuda_scaling_results_r64/` - Ratio 64 results

Each directory contains:
- `gpu_info.txt` - GPU hardware information
- `cuda_grid_*_mfem.txt` - Individual CUDA test results
- `omp_grid_*_mfem.txt` - Individual OpenMP test results
- `SUMMARY.txt` - Summary report with speedup calculations

## Test Matrix

Each script tests these mesh sizes:
- 16×16 (289 DOFs)
- 32×32 (1,089 DOFs)
- 64×64 (4,225 DOFs)
- 96×96 (9,409 DOFs)
- 128×128 (16,641 DOFs)
- 192×192 (37,249 DOFs)
- 256×256 (66,049 DOFs)
- 384×384 (148,225 DOFs)
- 512×512 (263,169 DOFs)

## Effective Fine Mesh Sizes

The effective fine mesh size = coarse mesh × ratio:

| Coarse Mesh | r=8 Fine | r=16 Fine | r=32 Fine | r=64 Fine |
|-------------|----------|-----------|-----------|-----------|
| 16×16       | 128×128  | 256×256   | 512×512   | 1024×1024 |
| 32×32       | 256×256  | 512×512   | 1024×1024 | 2048×2048 |
| 64×64       | 512×512  | 1024×1024 | 2048×2048 | 4096×4096 |
| 128×128     | 1024×1024| 2048×2048 | 4096×4096 | 8192×8192 |
| 256×256     | 2048×2048| 4096×4096 | 8192×8192 | 16384×16384 |
| 512×512     | 4096×4096| 8192×8192 | 16384×16384 | 32768×32768 |

## Example Workflow

```bash
# Complete workflow for all ratios

# Ratio 8
vi src/ScaledLaplacian.h  # Set ratio = 8
cd build && make cuda_assembly openmp_assembly -j8 && cd ..
./run_cpp_cuda_r8.sh

# Ratio 16
vi src/ScaledLaplacian.h  # Set ratio = 16
cd build && make cuda_assembly openmp_assembly -j8 && cd ..
./run_cpp_cuda_r16.sh

# Ratio 32
vi src/ScaledLaplacian.h  # Set ratio = 32
cd build && make cuda_assembly openmp_assembly -j8 && cd ..
./run_cpp_cuda_r32.sh

# Ratio 64
vi src/ScaledLaplacian.h  # Set ratio = 64
cd build && make cuda_assembly openmp_assembly -j8 && cd ..
./run_cpp_cuda_r64.sh
```

## Analyzing Results

After all runs complete, use the summary reports:

```bash
cat cpp_cuda_scaling_results_r8/SUMMARY.txt
cat cpp_cuda_scaling_results_r16/SUMMARY.txt
cat cpp_cuda_scaling_results_r32/SUMMARY.txt
cat cpp_cuda_scaling_results_r64/SUMMARY.txt
```

Compare:
- GPU speedup vs OpenMP for each ratio
- How speedup changes with mesh size
- Which ratio provides best GPU performance
- Memory usage and scalability limits

## Troubleshooting

**Script says "ratio mismatch":**
- You didn't recompile after changing the ratio
- Edit `src/ScaledLaplacian.h` and rebuild

**CUDA out of memory:**
- Reduce maximum mesh size (512×512 may be too large for some GPUs with ratio=64)
- Monitor GPU memory with `nvidia-smi` during runs

**OpenMP tests missing:**
- `openmp_assembly.exe` not built
- Rebuild with: `make openmp_assembly`

## Performance Tips

1. **Warm up GPU**: First run may be slower due to JIT compilation
2. **System load**: Run on dedicated GPU node for consistent results
3. **Thread count**: OpenMP uses 8 threads by default, adjust if needed
4. **Batch mode**: Run all tests in one session to minimize setup overhead

## Comparison with Julia

Julia CUDA results are in `julia_cuda_scaling_results/`. Compare:
- C++ vs Julia performance for same problem sizes
- Different ratio strategies
- Memory efficiency
- Code complexity vs performance trade-offs
