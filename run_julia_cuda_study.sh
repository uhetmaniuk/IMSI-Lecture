#!/bin/bash
#
# GPU (CUDA) Scaling Study for MFEM Assembly
# Tests problem size scaling and ratio scaling using mfem_cuda_block.jl
#

OUTPUT_DIR="julia_cuda_scaling_results"
mkdir -p $OUTPUT_DIR

echo "=========================================================================="
echo "MFEM Assembly GPU Scaling Study - Julia CUDA (Block Kernel)"
echo "=========================================================================="
echo ""
echo "This script studies:"
echo "  1. Grid scaling (varying mesh size, fixed ratio=32)"
echo "  2. Ratio scaling (varying ratio, fixed mesh=128×128)"
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

# ==============================================================================
# Study 1: Grid Scaling (varying mesh size, fixed ratio=32)
# ==============================================================================
echo "=========================================================================="
echo "Study 1: Grid Scaling (fixed ratio=32)"
echo "=========================================================================="
echo "Mesh sizes: 16×16, 32×32, 64×64, 128×128, 256×256, 512×512"
echo ""

RATIO=32
for nx in 16 32 64 128 256 512; do
    echo "  Running ${nx}×${nx} mesh with ratio=${RATIO}..."
    output_file="$OUTPUT_DIR/grid_${nx}x${nx}_r${RATIO}.txt"

    julia julia/mfem_cuda_block.jl $nx $nx $RATIO > $output_file 2>&1

    # Extract key metrics from mfem_cuda_block.jl output
    assembly_time=$(grep "MFEM assembly:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    solve_time=$(grep "Solve:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    reconstruct_time=$(grep "Reconstruction:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    total_time=$(grep "Total:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    dofs=$(grep "DOFs (coarse):" $output_file | awk '{print $NF}')
    fine_res=$((nx * RATIO))

    echo "    DOFs: ${dofs}, Fine: ${fine_res}×${fine_res}, Assembly: ${assembly_time} ms, Solve: ${solve_time} ms, Total: ${total_time} ms"
done

echo ""
echo "Grid scaling complete!"
echo ""

# ==============================================================================
# Study 2: Ratio Scaling (varying ratio, fixed mesh=128×128)
# ==============================================================================
echo "=========================================================================="
echo "Study 2: Ratio Scaling (fixed mesh=128×128)"
echo "=========================================================================="
echo "Ratios: 8, 16, 32"
echo ""

NX=128
NY=128
for ratio in 8 16 32; do
    echo "  Running ${NX}×${NY} mesh with ratio=${ratio}..."
    output_file="$OUTPUT_DIR/ratio_${ratio}_mesh${NX}x${NY}.txt"

    julia julia/mfem_cuda_block.jl $NX $NY $ratio > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "MFEM assembly:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    solve_time=$(grep "Solve:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    reconstruct_time=$(grep "Reconstruction:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    total_time=$(grep "Total:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    fine_res=$((NX * ratio))

    echo "    Fine res: ${fine_res}×${fine_res}, Assembly: ${assembly_time} ms, Solve: ${solve_time} ms, Total: ${total_time} ms"
done

echo ""
echo "Ratio scaling complete!"
echo ""

# ==============================================================================
# Generate Summary Report
# ==============================================================================
echo "=========================================================================="
echo "Generating Summary Report"
echo "=========================================================================="

summary_file="$OUTPUT_DIR/SUMMARY.txt"
{
    echo "=========================================================================="
    echo "MFEM GPU (CUDA) Scaling Study - Summary Report (mfem_cuda_block.jl)"
    echo "=========================================================================="
    echo ""
    echo "Generated: $(date)"
    echo ""

    echo "--- Grid Scaling (ratio=32) ---"
    echo "Mesh      | DOFs    | Fine Res      | Assembly (ms) | Solve (ms) | Total (ms)"
    echo "----------|---------|---------------|---------------|------------|------------"

    for nx in 16 32 64 128 256 512; do
        output_file="$OUTPUT_DIR/grid_${nx}x${nx}_r32.txt"
        if [ -f "$output_file" ]; then
            dofs=$(grep "DOFs (coarse):" $output_file | awk '{print $NF}')
            fine_res=$((nx * 32))
            assembly=$(grep "MFEM assembly:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            solve=$(grep "Solve:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            total=$(grep "Total:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            printf "%4dx%-5d | %7s | %5dx%-7d | %13s | %10s | %10s\n" \
                   $nx $nx "$dofs" $fine_res $fine_res "$assembly" "$solve" "$total"
        fi
    done
    echo ""

    echo "--- Ratio Scaling (128×128 mesh) ---"
    echo "Ratio | Fine Res      | Assembly (ms) | Solve (ms) | Total (ms)"
    echo "------|---------------|---------------|------------|------------"

    for ratio in 8 16 32; do
        output_file="$OUTPUT_DIR/ratio_${ratio}_mesh128x128.txt"
        if [ -f "$output_file" ]; then
            fine_res=$((128 * ratio))
            assembly=$(grep "MFEM assembly:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            solve=$(grep "Solve:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            total=$(grep "Total:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            printf "%5d | %5dx%-7d | %13s | %10s | %10s\n" \
                   $ratio $fine_res $fine_res "$assembly" "$solve" "$total"
        fi
    done
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
echo "=========================================================================="
