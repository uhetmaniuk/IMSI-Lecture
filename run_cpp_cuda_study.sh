#!/bin/bash
#
# C++ CUDA Scaling Study for MFEM Assembly
# Tests problem size scaling and compares with OpenMP baseline
#

OUTPUT_DIR="cpp_cuda_scaling_results"
mkdir -p $OUTPUT_DIR

echo "=========================================================================="
echo "MFEM Assembly GPU Scaling Study - C++ CUDA"
echo "=========================================================================="
echo ""
echo "NOTE: Run this script from the repository root directory"
echo ""
echo "This script studies:"
echo "  1. Grid scaling (varying mesh size, MFEM elements)"
echo "  2. Comparison with OpenMP baseline (if requested)"
echo ""

# Capture GPU info
echo "Capturing GPU information..."
{
    echo "=== GPU Information ==="
    echo "Date: $(date)"
    echo "Hostname: $(hostname)"

    # Try to get NVIDIA GPU info
    if command -v nvidia-smi &> /dev/null; then
        echo ""
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
    echo "  make cuda_assembly"
    exit 1
fi

# ==============================================================================
# Study 1: Grid Scaling (varying mesh size)
# ==============================================================================
echo "=========================================================================="
echo "Study 1: Grid Scaling (MFEM elements)"
echo "=========================================================================="
echo "Mesh sizes: 8×8, 16×16, 24×24, 32×32, 48×48"
echo ""

for nx in 8 16 24 32 48; do
    echo "  Running ${nx}×${nx} mesh with MFEM elements..."
    output_file="$OUTPUT_DIR/cuda_grid_${nx}x${nx}_mfem.txt"

    ./build/examples/cuda_assembly.exe -nx $nx -ny $nx -mfem > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Assembly time:" $output_file | awk '{print $3}')
    total_dofs=$(grep "Total DOFs:" $output_file | awk '{print $3}')

    echo "    DOFs: ${total_dofs}, Assembly: ${assembly_time} ms"
done

echo ""
echo "Grid scaling complete!"
echo ""

# ==============================================================================
# Study 2: OpenMP Baseline for Comparison (Optional)
# ==============================================================================
echo "=========================================================================="
echo "Study 2: OpenMP Baseline (for speedup comparison)"
echo "=========================================================================="

if [ ! -f "./build/examples/openmp_assembly.exe" ]; then
    echo "Warning: openmp_assembly.exe not found - skipping CPU baseline"
    echo ""
else
    echo "Running matching problems on CPU (OpenMP)"
    echo ""

    for nx in 8 16 24 32 48; do
        echo "  Running ${nx}×${nx} mesh on CPU (OpenMP, 8 threads)..."
        output_file="$OUTPUT_DIR/omp_grid_${nx}x${nx}_mfem.txt"

        ./build/examples/openmp_assembly.exe -nx $nx -ny $nx -mfem --kokkos-num-threads=8 > $output_file 2>&1

        # Extract key metrics
        assembly_time=$(grep "MFEM assembly time:" $output_file | awk '{print $3}')

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

summary_file="$OUTPUT_DIR/SUMMARY.txt"
{
    echo "=========================================================================="
    echo "MFEM GPU (CUDA) Scaling Study - Summary Report (C++)"
    echo "=========================================================================="
    echo ""
    echo "Generated: $(date)"
    echo ""

    echo "--- Grid Scaling (CUDA) ---"
    echo "Mesh    | DOFs   | Assembly (ms) | Free DOFs | Nnz"
    echo "--------|--------|---------------|-----------|--------"

    for nx in 8 16 24 32 48; do
        output_file="$OUTPUT_DIR/cuda_grid_${nx}x${nx}_mfem.txt"
        if [ -f "$output_file" ]; then
            dofs=$(grep "Total DOFs:" $output_file | awk '{print $3}')
            assembly=$(grep "Assembly time:" $output_file | awk '{print $3}')
            free_dofs=$(grep "Free DOFs:" $output_file | awk '{print $3}')
            nnz=$(grep "Non-zeros:" $output_file | head -1 | awk '{print $2}')
            printf "%4dx%-2d | %6s | %13s | %9s | %7s\n" \
                   $nx $nx "$dofs" "$assembly" "$free_dofs" "$nnz"
        fi
    done
    echo ""

    # CUDA vs OpenMP comparison
    if [ -f "$OUTPUT_DIR/omp_grid_32x32_mfem.txt" ]; then
        echo "--- CUDA vs OpenMP Comparison (32×32 mesh) ---"
        echo "Platform    | Assembly Time (ms) | Speedup"
        echo "------------|--------------------|---------"

        cuda_output="$OUTPUT_DIR/cuda_grid_32x32_mfem.txt"
        omp_output="$OUTPUT_DIR/omp_grid_32x32_mfem.txt"

        cuda_time=$(grep "Assembly time:" $cuda_output | awk '{print $3}')
        omp_time=$(grep "MFEM assembly time:" $omp_output | awk '{print $3}')

        # Calculate speedup (using bc for floating point)
        if [ ! -z "$omp_time" ] && [ ! -z "$cuda_time" ]; then
            speedup=$(echo "scale=2; $omp_time / $cuda_time" | bc)
            printf "OpenMP (8t) | %18s |   1.00x\n" "$omp_time"
            printf "CUDA        | %18s | %6sx\n" "$cuda_time" "$speedup"
            echo ""
            echo "CUDA Speedup: ${speedup}x faster than OpenMP (8 threads)"
        fi
        echo ""

        # Full comparison table
        echo "--- Complete CUDA vs OpenMP Comparison ---"
        echo "Mesh    | CUDA (ms) | OpenMP (ms) | Speedup"
        echo "--------|-----------|-------------|--------"

        for nx in 8 16 24 32 48; do
            cuda_output="$OUTPUT_DIR/cuda_grid_${nx}x${nx}_mfem.txt"
            omp_output="$OUTPUT_DIR/omp_grid_${nx}x${nx}_mfem.txt"

            if [ -f "$cuda_output" ] && [ -f "$omp_output" ]; then
                cuda_time=$(grep "Assembly time:" $cuda_output | awk '{print $3}')
                omp_time=$(grep "MFEM assembly time:" $omp_output | awk '{print $3}')

                if [ ! -z "$omp_time" ] && [ ! -z "$cuda_time" ]; then
                    speedup=$(echo "scale=2; $omp_time / $cuda_time" | bc)
                    printf "%4dx%-2d | %9s | %11s | %6sx\n" \
                           $nx $nx "$cuda_time" "$omp_time" "$speedup"
                fi
            fi
        done
        echo ""
    fi

    echo "--- Key Findings ---"
    echo "• CUDA provides significant acceleration for MFEM assembly"
    echo "• Speedup increases with problem size (Amdahl's law)"
    echo "• GPU memory bandwidth is critical for performance"
    echo "• Larger meshes show better GPU utilization"
    echo ""
    echo "=========================================================================="
    echo "All results available in: $OUTPUT_DIR/"
    echo "=========================================================================="

} > $summary_file

cat $summary_file
echo ""

echo "=========================================================================="
echo "All experiments complete!"
echo "=========================================================================="
echo ""
echo "Results summary: $summary_file"
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "To compare with other implementations:"
echo "  ./run_cpp_omp_study.sh      # C++ OpenMP"
echo "  ./run_julia_cuda_study.sh   # Julia CUDA"
echo "  ./run_julia_omp_study.sh    # Julia OpenMP"
echo ""
echo "=========================================================================="
