#!/bin/bash
#
# CPU (OpenMP) Scaling Study for MFEM Assembly (Julia)
# Tests thread scaling, grid scaling, and ratio scaling
#

# Configuration
RATIO=8

OUTPUT_DIR="julia_omp_scaling_results_r${RATIO}"
mkdir -p $OUTPUT_DIR

echo "=========================================================================="
echo "MFEM Assembly CPU Scaling Study - Julia OpenMP"
echo "=========================================================================="
echo ""
echo "This script studies:"
echo "  1. Thread scaling (varying thread count, fixed problem)"
echo "  2. Grid scaling (varying mesh size, fixed ratio)"
echo "  3. Ratio scaling (varying ratio, fixed mesh)"
echo ""

# Capture system info
echo "Capturing system information..."
{
    echo "=== System Information ==="
    echo "Date: $(date)"
    echo "Hostname: $(hostname)"
    echo "Julia version:"
    julia --version
    echo ""
} > $OUTPUT_DIR/system_info.txt

cat $OUTPUT_DIR/system_info.txt
echo ""

# ==============================================================================
# Study 1: Thread Scaling (fixed problem size, optimal configuration)
# ==============================================================================
echo "=========================================================================="
echo "Study 0: Thread Scaling (48x48 mesh, ratio=${RATIO})"
echo "=========================================================================="
echo "Thread counts: 1, 2, 4, 8, 12, 16"
echo ""

for threads in 1 2 4 8 12 16; do
    echo "  Running with $threads thread(s)..."
    output_file="$OUTPUT_DIR/thread_${threads}_48x48_r${RATIO}.txt"

    julia -t $threads julia/mfem_assembly.jl 48 48 $RATIO > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    echo "    Assembly: ${assembly_time} ms"
done

echo "=========================================================================="
echo "Study 1: Thread Scaling (128×128 mesh, ratio=${RATIO})"
echo "=========================================================================="
echo "Thread counts: 1, 2, 4, 8, 12, 16"
echo ""

for threads in 1 2 4 8 12 16; do
    echo "  Running with $threads thread(s)..."
    output_file="$OUTPUT_DIR/thread_${threads}_128x128_r${RATIO}.txt"

    julia -t $threads julia/mfem_assembly.jl 128 128 $RATIO > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    echo "    Assembly: ${assembly_time} ms"
done

echo ""
echo "Thread scaling complete!"
echo ""

# ==============================================================================
# Study 2: Grid Scaling (varying mesh size, fixed ratio, fixed threads=8)
# ==============================================================================
echo "=========================================================================="
echo "Study 2: Grid Scaling (ratio=${RATIO}, 8 threads)"
echo "=========================================================================="
echo "Mesh sizes: 16×16, 32×32, 48x48, 64×64, 128×128, 256×256"
echo ""

THREADS=8
for nx in 16 32 48 64 128 256; do
    echo "  Running ${nx}×${nx} mesh..."
    output_file="$OUTPUT_DIR/grid_${nx}x${nx}_r${RATIO}_t${THREADS}.txt"

    julia -t $THREADS julia/mfem_assembly.jl $nx $nx $RATIO > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    dofs=$(grep "DOFs (coarse):" $output_file | awk '{print $NF}')

    echo "    DOFs: ${dofs}, Assembly: ${assembly_time} ms"
done

echo ""
echo "Grid scaling complete!"
echo ""

# ==============================================================================
# Study 3: Ratio Scaling (varying ratio, fixed mesh=128×128, fixed threads=8)
# ==============================================================================
echo "=========================================================================="
echo "Study 3: Ratio Scaling (128×128 mesh, 8 threads)"
echo "=========================================================================="
echo "Ratios: 4, 8, 16, 32, 64"
echo ""

NX=128
NY=128
THREADS=8
for ratio in 4 8 16 32 64; do
    echo "  Running ratio=${ratio}..."
    output_file="$OUTPUT_DIR/ratio_${ratio}_mesh${NX}x${NY}_t${THREADS}.txt"

    julia -t $THREADS julia/mfem_assembly.jl $NX $NY $ratio > $output_file 2>&1

    # Extract key metrics
    assembly_time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
    fine_res=$((NX * ratio))

    echo "    Fine res: ${fine_res}×${fine_res}, Assembly: ${assembly_time} ms"
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
    echo "MFEM CPU (OpenMP) Scaling Study - Summary Report (Julia)"
    echo "=========================================================================="
    echo ""
    echo "Generated: $(date)"
    echo ""

    echo "--- Thread Scaling (48x48 mesh, ratio=${RATIO}) ---"
    echo "Threads | Assembly (ms) | Speedup | Efficiency"
    echo "--------|---------------|---------|------------"

    # Get baseline (1 thread) time
    baseline_file="$OUTPUT_DIR/thread_1_48x48_r${RATIO}.txt"
    if [ -f "$baseline_file" ]; then
        baseline_time=$(grep "Total assembly time:" $baseline_file | grep -oE '[0-9]+\.[0-9]+' | head -1)

        for threads in 1 2 4 8 12 16; do
            output_file="$OUTPUT_DIR/thread_${threads}_48x48_r${RATIO}.txt"
            if [ -f "$output_file" ]; then
                time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)

                if [ ! -z "$time" ] && [ ! -z "$baseline_time" ]; then
                    speedup=$(echo "scale=2; $baseline_time / $time" | bc)
                    efficiency=$(echo "scale=1; 100 * $speedup / $threads" | bc)
                    printf "%7d | %13s | %7sx | %10s%%\n" \
                           $threads "$time" "$speedup" "$efficiency"
                fi
            fi
        done
    fi
    echo ""

    echo "--- Thread Scaling (128×128 mesh, ratio=${RATIO}) ---"
    echo "Threads | Assembly (ms) | Speedup | Efficiency"
    echo "--------|---------------|---------|------------"

    # Get baseline (1 thread) time
    baseline_file="$OUTPUT_DIR/thread_1_128x128_r${RATIO}.txt"
    if [ -f "$baseline_file" ]; then
        baseline_time=$(grep "Total assembly time:" $baseline_file | grep -oE '[0-9]+\.[0-9]+' | head -1)

        for threads in 1 2 4 8 12 16; do
            output_file="$OUTPUT_DIR/thread_${threads}_128x128_r${RATIO}.txt"
            if [ -f "$output_file" ]; then
                time=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)

                if [ ! -z "$time" ] && [ ! -z "$baseline_time" ]; then
                    speedup=$(echo "scale=2; $baseline_time / $time" | bc)
                    efficiency=$(echo "scale=1; 100 * $speedup / $threads" | bc)
                    printf "%7d | %13s | %7sx | %10s%%\n" \
                           $threads "$time" "$speedup" "$efficiency"
                fi
            fi
        done
    fi
    echo ""

    echo "--- Grid Scaling (ratio=${RATIO}, 8 threads) ---"
    echo "Mesh    | DOFs   | Fine Res   | Assembly (ms)"
    echo "--------|--------|------------|---------------"

    for nx in 16 32 48 64 128 256; do
        output_file="$OUTPUT_DIR/grid_${nx}x${nx}_r${RATIO}_t8.txt"
        if [ -f "$output_file" ]; then
            dofs=$(grep "DOFs (coarse):" $output_file | awk '{print $NF}')
            fine_res=$((nx * RATIO))
            assembly=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            printf "%4dx%-3d | %6s | %4dx%-5d | %13s\n" \
                   $nx $nx "$dofs" $fine_res $fine_res "$assembly"
        fi
    done
    echo ""

    echo "--- Ratio Scaling (128×128 mesh, 8 threads) ---"
    echo "Ratio | Fine Res    | Assembly (ms)"
    echo "------|-------------|---------------"

    for ratio in 4 8 16 32 64; do
        output_file="$OUTPUT_DIR/ratio_${ratio}_mesh128x128_t8.txt"
        if [ -f "$output_file" ]; then
            fine_res=$((128 * ratio))
            assembly=$(grep "Total assembly time:" $output_file | grep -oE '[0-9]+\.[0-9]+' | head -1)
            printf "%5d | %5dx%-6d | %13s\n" \
                   $ratio $fine_res $fine_res "$assembly"
        fi
    done
    echo ""

    echo "--- Key Findings ---"
    echo "• OpenMP scales efficiently with thread count (check efficiency %)"
    echo "• Assembly time scales with problem size (check grid scaling)"
    echo "• Higher ratios increase local problem size (check ratio scaling)"
    echo "• Compare with GPU results for hardware acceleration insights"
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
echo "To compare with GPU implementation:"
echo "  ./run_julia_cuda_study.sh"
echo ""
echo "=========================================================================="
