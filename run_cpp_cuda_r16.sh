#!/bin/bash
#
# C++ CUDA Scaling Study for MFEM Assembly - RATIO=16
# Tests problem size scaling and compares with OpenMP baseline
#
# IMPORTANT: Before running this script, you must recompile with ratio=16
#   Edit src/ScaledLaplacian.h line 144: static constexpr int ratio = 16;
#   Then: cd build && cmake .. && make cuda_assembly openmp_assembly
#

OUTPUT_DIR="cpp_cuda_scaling_results_r16"
RATIO=16

echo "=========================================================================="
echo "MFEM Assembly GPU Scaling Study - C++ CUDA (RATIO=$RATIO)"
echo "=========================================================================="
echo ""
echo "IMPORTANT: This script assumes ratio=$RATIO in src/ScaledLaplacian.h"
echo "           If you haven't recompiled with ratio=$RATIO, results will be incorrect!"
echo ""
read -p "Have you recompiled with ratio=$RATIO? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Please recompile first:"
    echo "  1. Edit src/ScaledLaplacian.h line 144: static constexpr int ratio = $RATIO;"
    echo "  2. cd build && cmake .. && make cuda_assembly openmp_assembly"
    exit 1
fi
echo ""

mkdir -p $OUTPUT_DIR

echo "NOTE: Run this script from the repository root directory"
echo ""
echo "This script studies:"
echo "  1. Grid scaling (varying mesh size, MFEM elements)"
echo "  2. Comparison with OpenMP baseline"
echo "  3. MFEM ratio: $RATIO (fine mesh is ${RATIO}x${RATIO} per coarse element)"
echo ""

# Capture GPU info
echo "Capturing GPU information..."
{
    echo "=== GPU Information ==="
    echo "Date: $(date)"
    echo "Hostname: $(hostname)"
    echo "MFEM Ratio: $RATIO"
    echo ""

    # Try to get NVIDIA GPU info
    if command -v nvidia-smi &> /dev/null; then
        nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
    else
        echo "nvidia-smi not found - unable to query GPU details"
    fi
    echo ""
} > $OUTPUT_DIR/gpu_info.txt

cat $OUTPUT_DIR/gpu_info.txt
echo ""

# Check if executables exist
if [ ! -f "./build/examples/cuda_assembly.exe" ]; then
    echo "Error: cuda_assembly.exe not found in build/examples/"
    echo "Please build with:"
    echo "  cd build"
    echo "  cmake .. -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_<YOUR_GPU>=ON"
    echo "  make cuda_assembly openmp_assembly"
    exit 1
fi

# ==============================================================================
# Study 1: Grid Scaling (varying mesh size) - CUDA
# ==============================================================================
echo "=========================================================================="
echo "Study 1: Grid Scaling - CUDA (MFEM elements, ratio=$RATIO)"
echo "=========================================================================="
echo "Mesh sizes: 16×16, 32×32, 64×64, 96×96, 128×128, 192×192, 256×256, 384×384, 512×512"
echo ""

for nx in 16 32 64 96 128 192 256 384 512; do
    echo "  Running ${nx}×${nx} mesh on CUDA..."
    output_file="$OUTPUT_DIR/cuda_grid_${nx}x${nx}_mfem.txt"

    ./build/examples/cuda_assembly.exe -nx $nx -ny $nx -mfem > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    total_dofs=$(grep "Total DOFs:" $output_file | awk '{print $NF}')

    echo "    DOFs: ${total_dofs}, Assembly: ${assembly_time} ms"
done

echo ""
echo "CUDA grid scaling complete!"
echo ""

# ==============================================================================
# Study 2: OpenMP Baseline for Comparison
# ==============================================================================
echo "=========================================================================="
echo "Study 2: OpenMP Baseline (for speedup comparison)"
echo "=========================================================================="

if [ ! -f "./build/examples/openmp_assembly.exe" ]; then
    echo "Warning: openmp_assembly.exe not found - skipping CPU baseline"
    echo ""
else
    echo "Running matching problems on CPU (OpenMP, 8 threads)"
    echo ""

    for nx in 16 32 64 96 128 192 256 384 512; do
        echo "  Running ${nx}×${nx} mesh on OpenMP (8 threads)..."
        output_file="$OUTPUT_DIR/omp_grid_${nx}x${nx}_mfem.txt"

        ./build/examples/openmp_assembly.exe -nx $nx -ny $nx -mfem --kokkos-num-threads=8 > $output_file 2>&1

        assembly_time=$(grep "MFEM assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)

        echo "    Assembly: ${assembly_time} ms"
    done

    echo ""
    echo "OpenMP baseline complete!"
    echo ""
fi

# ==============================================================================
# Generate Summary Report
# ==============================================================================
echo "=========================================================================="
echo "Generating Summary Report"
echo "=========================================================================="

# Use Python script to generate summary
if [ -f "fix_cpp_summaries.py" ]; then
    python3 - <<EOF
from pathlib import Path
import sys
sys.path.insert(0, '.')
from fix_cpp_summaries import generate_summary
generate_summary('$OUTPUT_DIR', $RATIO)
EOF
    cat $OUTPUT_DIR/SUMMARY.txt
else
    echo "Warning: fix_cpp_summaries.py not found - cannot generate summary"
fi

echo ""
echo "=========================================================================="
echo "All experiments complete!"
echo "=========================================================================="
echo ""
echo "Results summary: $OUTPUT_DIR/SUMMARY.txt"
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "To run other ratios:"
echo "  ./run_cpp_cuda_r8.sh     # ratio=8"
echo "  ./run_cpp_cuda_r32.sh    # ratio=32"
echo "  ./run_cpp_cuda_r64.sh    # ratio=64"
echo ""
echo "=========================================================================="
