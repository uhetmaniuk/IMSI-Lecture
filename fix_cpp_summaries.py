#!/usr/bin/env python3
"""
Parse C++ CUDA scaling results and regenerate SUMMARY.txt files
"""

import re
import os
from pathlib import Path
from datetime import datetime

def parse_timing_file(filepath):
    """Extract timing information from a result file"""
    data = {}
    with open(filepath, 'r') as f:
        content = f.read()

        # Extract mesh size
        mesh_match = re.search(r'Generating mesh: (\d+) x (\d+)', content)
        if mesh_match:
            data['mesh'] = f"{mesh_match.group(1)}x{mesh_match.group(2)}"

        # Extract DOFs
        dof_match = re.search(r'Number of nodes:\s+(\d+)', content)
        if dof_match:
            data['dofs'] = int(dof_match.group(1))

        # Extract free DOFs
        free_dof_match = re.search(r'Free DOFs:\s+(\d+)', content)
        if free_dof_match:
            data['free_dofs'] = int(free_dof_match.group(1))

        # Extract nnz
        nnz_match = re.search(r'Non-zeros:\s+(\d+)', content)
        if nnz_match:
            data['nnz'] = int(nnz_match.group(1))

        # Extract ratio
        ratio_match = re.search(r'ratio=(\d+)', content)
        if ratio_match:
            data['ratio'] = int(ratio_match.group(1))

        # Extract assembly time (try both formats)
        asm_match = re.search(r'MFEM assembly time:\s+([\d.]+)\s+ms', content)
        if not asm_match:
            asm_match = re.search(r'^Assembly time:\s+([\d.]+)\s+ms', content, re.MULTILINE)
        if asm_match:
            data['assembly_time'] = float(asm_match.group(1))

    return data

def generate_summary(results_dir, ratio):
    """Generate SUMMARY.txt for a given results directory"""

    results_path = Path(results_dir)

    # Parse all CUDA files
    cuda_results = {}
    for cuda_file in sorted(results_path.glob('cuda_grid_*_mfem.txt')):
        data = parse_timing_file(cuda_file)
        if 'mesh' in data:
            cuda_results[data['mesh']] = data

    # Parse all OpenMP files
    omp_results = {}
    for omp_file in sorted(results_path.glob('omp_grid_*_mfem.txt')):
        data = parse_timing_file(omp_file)
        if 'mesh' in data:
            omp_results[data['mesh']] = data

    # Read GPU info
    gpu_info = "Unknown GPU"
    gpu_info_file = results_path / 'gpu_info.txt'
    if gpu_info_file.exists():
        with open(gpu_info_file, 'r') as f:
            for line in f:
                # Try NVIDIA format: "NVIDIA H100 PCIe, driver, memory"
                if 'NVIDIA' in line and ',' in line:
                    parts = line.split(',')
                    gpu_info = parts[0].strip()
                    break
                # Try older format
                if 'Device 0:' in line:
                    gpu_info = line.split('Device 0:')[1].strip()
                    break

    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 74)
    summary_lines.append("MFEM GPU (CUDA) Scaling Study - Summary Report (C++)")
    summary_lines.append("=" * 74)
    summary_lines.append("")
    summary_lines.append(f"Generated: {datetime.utcnow().strftime('%a %b %d %H:%M:%S UTC %Y')}")
    summary_lines.append(f"GPU: {gpu_info}")
    summary_lines.append(f"MFEM Ratio: {ratio}")
    summary_lines.append("")

    # Grid scaling table - CUDA
    summary_lines.append("--- Grid Scaling (CUDA) ---")
    summary_lines.append("Mesh    | DOFs   | Assembly (ms) | Free DOFs | Nnz")
    summary_lines.append("--------|--------|---------------|-----------|--------")

    for mesh in sorted(cuda_results.keys(), key=lambda x: int(x.split('x')[0])):
        data = cuda_results[mesh]
        summary_lines.append(
            f"{mesh:>7} | {data['dofs']:>6} | "
            f"{data['assembly_time']:>13.3f} | "
            f"{data.get('free_dofs', 0):>9} | "
            f"{data.get('nnz', 0):>7}"
        )
    summary_lines.append("")

    # Grid scaling table - OpenMP
    summary_lines.append("--- Grid Scaling (OpenMP, 8 threads) ---")
    summary_lines.append("Mesh    | DOFs   | Assembly (ms) | Free DOFs | Nnz")
    summary_lines.append("--------|--------|---------------|-----------|--------")

    for mesh in sorted(omp_results.keys(), key=lambda x: int(x.split('x')[0])):
        data = omp_results[mesh]
        summary_lines.append(
            f"{mesh:>7} | {data['dofs']:>6} | "
            f"{data['assembly_time']:>13.2f} | "
            f"{data.get('free_dofs', 0):>9} | "
            f"{data.get('nnz', 0):>7}"
        )
    summary_lines.append("")

    # Speedup comparison
    summary_lines.append("--- CUDA vs OpenMP Speedup ---")
    summary_lines.append("Mesh    | CUDA (ms) | OpenMP (ms) | Speedup")
    summary_lines.append("--------|-----------|-------------|--------")

    for mesh in sorted(cuda_results.keys(), key=lambda x: int(x.split('x')[0])):
        if mesh in omp_results:
            cuda_time = cuda_results[mesh]['assembly_time']
            omp_time = omp_results[mesh]['assembly_time']
            speedup = omp_time / cuda_time
            summary_lines.append(
                f"{mesh:>7} | {cuda_time:>9.2f} | {omp_time:>11.2f} | {speedup:>6.2f}x"
            )
    summary_lines.append("")

    # Pick a reference mesh for detailed comparison
    ref_mesh = '32x32'
    if ref_mesh in cuda_results and ref_mesh in omp_results:
        cuda_time = cuda_results[ref_mesh]['assembly_time']
        omp_time = omp_results[ref_mesh]['assembly_time']
        speedup = omp_time / cuda_time

        summary_lines.append(f"--- Detailed Comparison ({ref_mesh} mesh) ---")
        summary_lines.append(f"Platform    | Assembly Time (ms) | Speedup")
        summary_lines.append(f"------------|--------------------|---------")
        summary_lines.append(f"OpenMP (8t) | {omp_time:>18.2f} |   1.00x")
        summary_lines.append(f"CUDA        | {cuda_time:>18.2f} | {speedup:>6.2f}x")
        summary_lines.append("")

        if speedup > 1.0:
            summary_lines.append(f"CUDA Speedup: {speedup:.2f}x faster than OpenMP (8 threads)")
        else:
            slowdown = 1.0 / speedup
            summary_lines.append(f"CUDA Performance: {slowdown:.2f}x slower than OpenMP (8 threads)")

    summary_lines.append("")
    summary_lines.append("--- Key Findings ---")

    # Calculate average speedup
    speedups = []
    for mesh in cuda_results.keys():
        if mesh in omp_results:
            speedup = omp_results[mesh]['assembly_time'] / cuda_results[mesh]['assembly_time']
            speedups.append(speedup)

    if speedups:
        avg_speedup = sum(speedups) / len(speedups)
        if avg_speedup > 1.0:
            summary_lines.append(f"• Average CUDA speedup: {avg_speedup:.2f}x")
            summary_lines.append(f"• CUDA provides significant acceleration for MFEM assembly")
        else:
            summary_lines.append(f"• Average CUDA performance: {1.0/avg_speedup:.2f}x slower than OpenMP")
            summary_lines.append(f"• Small problem sizes may not benefit from GPU acceleration")

    summary_lines.append(f"• Problem sizes tested: {min(cuda_results.keys(), key=lambda x: int(x.split('x')[0]))} to {max(cuda_results.keys(), key=lambda x: int(x.split('x')[0]))}")
    summary_lines.append(f"• MFEM ratio: {ratio} (fine mesh is {ratio}x{ratio} per coarse element)")
    summary_lines.append(f"• GPU memory bandwidth and kernel launch overhead affect small problems")

    summary_lines.append("")
    summary_lines.append("=" * 74)
    summary_lines.append(f"All results available in: {results_path.name}/")
    summary_lines.append("=" * 74)
    summary_lines.append("")

    # Write summary file
    summary_file = results_path / 'SUMMARY.txt'
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))

    print(f"Generated: {summary_file}")
    return summary_lines

if __name__ == '__main__':
    # Generate summaries for both directories
    print("Generating fixed SUMMARY.txt files...\n")

    summary1 = generate_summary('cpp_cuda_scaling_results', ratio=32)
    print()
    summary2 = generate_summary('cpp_cuda_scaling_results_64', ratio=64)

    print("\nDone! Summary files have been regenerated with correct timing data.")
