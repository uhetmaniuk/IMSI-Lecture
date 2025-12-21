#!/bin/bash
#
# Scaling study for MFEM assembly
# Tests thread scaling and grid refinement
#

OUTPUT_DIR="scaling_results"
mkdir -p $OUTPUT_DIR

echo "=========================================="
echo "MFEM Assembly Scaling Study"
echo "=========================================="
echo ""

# ============================================================
# Study 1: Thread Scaling (fixed problem size)
# ============================================================
echo "Study 1: Thread Scaling (nx=32, ny=32, ratio=32)"
echo "Running with 1, 2, 4 threads..."
echo ""

for threads in 1 2 4; do
    echo "  Running with $threads thread(s)..."
    julia -t $threads julia/mfem_assembly.jl 32 32 32 > $OUTPUT_DIR/thread_${threads}.txt 2>&1
done

echo "Thread scaling study complete."
echo ""

# ============================================================
# Study 2: Grid Refinement (fixed threads)
# ============================================================
echo "Study 2: Grid Refinement (threads=2, ratio=32)"
echo "Running with mesh sizes: 16x16, 32x32, 48x48..."
echo ""

for nx in 16 32 48; do
    echo "  Running with ${nx}x${nx} mesh..."
    julia -t 2 julia/mfem_assembly.jl $nx $nx 32 > $OUTPUT_DIR/grid_${nx}.txt 2>&1
done

echo "Grid refinement study complete."
echo ""

# ============================================================
# Study 3: Ratio Scaling (fixed threads and coarse mesh)
# ============================================================
echo "Study 3: Ratio Scaling (threads=2, nx=32, ny=32)"
echo "Running with ratios: 16, 32, 64..."
echo ""

for ratio in 16 32 64; do
    echo "  Running with ratio=$ratio..."
    julia -t 2 julia/mfem_assembly.jl 32 32 $ratio > $OUTPUT_DIR/ratio_${ratio}.txt 2>&1
done

echo "Ratio scaling study complete."
echo ""

echo "=========================================="
echo "All experiments complete!"
echo "Results saved in: $OUTPUT_DIR/"
echo "=========================================="
