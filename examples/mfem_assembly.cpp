///
/// @file mfem_assembly.cpp
/// @brief Example demonstrating OpenMP-accelerated MFEM assembly using ScaledLaplacian
///
/// This example shows how to use the ScaledLaplacian class with Kokkos::OpenMP
/// execution space to assemble the Multiscale Finite Element Method (MFEM) operator.
///
/// Solves: -∇·(α∇u) = f on a 2D rectangular domain using MFEM_L elements
/// with homogeneous Dirichlet boundary conditions
///
/// MFEM uses static condensation: solves coarse problem, then reconstructs fine solution
///

#include <Kokkos_Core.hpp>
#include <cmath>
#include <functional>
#include <iostream>

#include "../src/Element.h"
#include "../src/Mesh.h"
#include "../src/MeshUtils.h"
#include "../src/ScaledLaplacian.h"
#include "../src/SymmetricSparse.hpp"
#include "../src/Utils.h"

using namespace IMSI;

// Define material coefficients and forcing term as functors at namespace scope
// These must be KOKKOS_INLINE_FUNCTION compatible for device execution
struct AxFunctor {
  KOKKOS_INLINE_FUNCTION
  double operator()(double x, double y, double z) const { return 1.0; }
};

struct AyFunctor {
  KOKKOS_INLINE_FUNCTION
  double operator()(double x, double y, double z) const { return 1.0; }
};

struct FFunctor {
  KOKKOS_INLINE_FUNCTION
  double operator()(double x, double y, double z) const {
    // Manufactured solution: u = sin(pi*x) * sin(pi*y)
    // => -Laplacian(u) = 2*pi^2 * sin(pi*x) * sin(pi*y)
    constexpr double pi = 3.14159265358979323846;
    using Kokkos::sin;  // Device-safe sin function
    return 2.0 * pi * pi * sin(pi * x) * sin(pi * y);
  }
};

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::cout << "=== OpenMP MFEM Assembly Example (using ScaledLaplacian) ===" << std::endl;
    std::cout << "Kokkos execution space: " << typeid(Kokkos::OpenMP).name() << std::endl;
    std::cout << std::endl;

    // Check if OpenMP is available
#ifdef KOKKOS_ENABLE_OPENMP
    std::cout << "OpenMP is enabled" << std::endl;
    std::cout << "OpenMP thread count: " << Kokkos::OpenMP().concurrency() << " threads" << std::endl;
#else
    std::cerr << "ERROR: This example requires Kokkos with OpenMP support" << std::endl;
    std::cerr << "       Rebuild with -DKokkos_ENABLE_OPENMP=ON" << std::endl;
    Kokkos::finalize();
    return 1;
#endif

    // ========================================================================
    // Problem setup - MFEM Configuration
    // ========================================================================

    const int coarseSize = 32;  // 32x32 coarse mesh
    const int ratio = 32;        // Fine mesh has 32x finer resolution

    // ========================================================================
    // Mesh generation
    // ========================================================================

    IMSI::DomainParams dParams;
    dParams.numElePerDir[0] = coarseSize;
    dParams.numElePerDir[1] = coarseSize;
    dParams.omega           = IMSI::DomainType::Rectangle;
    dParams.cellType        = IMSI::ElementType::MFEM_L;

    std::cout << "Generating MFEM mesh:" << std::endl;
    std::cout << "  Coarse mesh: " << dParams.numElePerDir[0] << " x " << dParams.numElePerDir[1]
              << " MFEM_L elements" << std::endl;
    std::cout << "  Ratio: " << ratio << std::endl;
    std::cout << "  Effective fine mesh: " << coarseSize * ratio << " x " << coarseSize * ratio
              << " (" << (coarseSize * ratio + 1) * (coarseSize * ratio + 1) << " fine nodes)" << std::endl;

    auto mesh = IMSI::GenerateMesh(dParams);
    std::cout << "  Number of coarse elements: " << mesh.NumberCells() << std::endl;
    std::cout << "  Number of coarse nodes:    " << mesh.NumberVertices() << std::endl;

    // ========================================================================
    // Build mesh connectivity and graph coloring
    // ========================================================================

    std::cout << "\nBuilding mesh connectivity..." << std::endl;
    auto meshConn = GetMeshConnectivity(mesh, true);  // true = use coloring

    std::cout << "  Number of colors: " << meshConn.c2e.numRows() << std::endl;
    for (int ic = 0; ic < meshConn.c2e.numRows(); ++ic) {
      std::cout << "    Color " << ic << ": " << meshConn.c2e.row_map(ic + 1) - meshConn.c2e.row_map(ic)
                << " elements" << std::endl;
    }

    // ========================================================================
    // Create sparse matrix structure
    // ========================================================================

    std::cout << "\nAllocating sparse matrix structure..." << std::endl;

    auto const numNodes = mesh.NumberVertices();
    auto const numNonZeros = meshConn.n2n.entries.extent(0);

    std::cout << "  Matrix size (coarse): " << numNodes << " x " << numNodes << std::endl;
    std::cout << "  Non-zeros:            " << numNonZeros << std::endl;

    // OpenMP execution space views (host-accessible)
    Kokkos::View<double*, Kokkos::OpenMP> rhs("rhs", numNodes);
    Kokkos::View<size_t*, Kokkos::OpenMP> matRowPtr("matRowPtr", numNodes + 1);
    Kokkos::View<int*, Kokkos::OpenMP> matColIdx("matColIdx", numNonZeros);
    Kokkos::View<double*, Kokkos::OpenMP> matValues("matValues", numNonZeros);

    // Copy graph structure from n2n (on host)
    for (size_t i = 0; i <= numNodes; ++i) {
      matRowPtr(i) = meshConn.n2n.row_map(i);
    }
    for (size_t i = 0; i < numNonZeros; ++i) {
      matColIdx(i) = meshConn.n2n.entries(i);
      matValues(i) = 0.0;
    }

    // Initialize to zero
    Kokkos::deep_copy(rhs, 0.0);

    // ========================================================================
    // Assembly with OpenMP - MFEM
    // ========================================================================

    std::cout << "\n=== Starting OpenMP MFEM Assembly ===" << std::endl;

    // Create coefficient functor instances
    AxFunctor ax_func;
    AyFunctor ay_func;
    FFunctor f_func;

    // Instantiate ScaledLaplacian with OpenMP ExecutionSpace and functor types
    auto scalarLap = ScaledLaplacian<Kokkos::OpenMP, AxFunctor, AyFunctor, FFunctor>(
        meshConn, RuleType::Gauss, 2, ax_func, ay_func, f_func);

    Kokkos::Timer timer;
    scalarLap.GetLinearSystemMFEM(rhs, matRowPtr, matColIdx, matValues);
    Kokkos::fence();
    double assemblyTime = timer.seconds();

    std::cout << "=== MFEM Assembly Complete ===" << std::endl;
    std::cout << "Assembly time: " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // Apply boundary conditions
    // ========================================================================

    std::cout << "Applying boundary conditions..." << std::endl;

    auto const&      bdyNodes = mesh.GetBoundaryNodes();
    std::vector<int> globalToFree(mesh.NumberVertices());
    std::vector<int> freeToGlobal(globalToFree.size() - bdyNodes.size());
    IMSI::MapDegreesOfFreedom(bdyNodes, globalToFree, freeToGlobal);

    auto numFreeDofs = freeToGlobal.size();
    std::cout << "  Total coarse DOFs:    " << numNodes << std::endl;
    std::cout << "  Boundary DOFs:        " << numNodes - numFreeDofs << std::endl;
    std::cout << "  Free coarse DOFs:     " << numFreeDofs << std::endl;

    // Build reduced system
    std::vector<int> newRowPtr(numFreeDofs + 1, 0);
    std::vector<int> newColIdx;
    std::vector<double> newValues;
    std::vector<double> newRhs(numFreeDofs);

    newRowPtr[0] = 0;
    for (size_t i = 0; i < numFreeDofs; ++i) {
      auto gDof = freeToGlobal[i];
      newRhs[i] = rhs(gDof);
      for (size_t k = matRowPtr(gDof); k < matRowPtr(gDof + 1); ++k) {
        auto gCol = matColIdx(k);
        if (globalToFree[gCol] != -1) {
          newColIdx.push_back(globalToFree[gCol]);
          newValues.push_back(matValues(k));
        }
      }
      newRowPtr[i + 1] = newColIdx.size();
    }

    std::cout << "  Reduced matrix size: " << numFreeDofs << " x " << numFreeDofs << std::endl;
    std::cout << "  Reduced non-zeros:   " << newColIdx.size() << std::endl;

    // ========================================================================
    // Solve on host (using CPU sparse solver)
    // ========================================================================

    std::cout << "\nSolving linear system on host..." << std::endl;

    SymmetricSparse<double> solver(
        numFreeDofs, newColIdx.size(), newRowPtr.data(), newColIdx.data(), newValues.data(), true);

    timer.reset();
    solver.factor();
    double factorTime = timer.seconds();
    std::cout << "  Factorization time: " << factorTime * 1000.0 << " ms" << std::endl;

    std::vector<double> solution(numFreeDofs);
    timer.reset();
    solver.Solve(1, newRhs.data(), solution.data());
    double solveTime = timer.seconds();
    std::cout << "  Solve time:         " << solveTime * 1000.0 << " ms" << std::endl;

    // ========================================================================
    // Output results - Coarse solution
    // ========================================================================

    std::cout << "\nWriting coarse solution to file..." << std::endl;

    std::vector<double> fullSolution(numNodes, 0.0);
    for (size_t i = 0; i < numFreeDofs; ++i) {
      fullSolution[freeToGlobal[i]] = solution[i];
    }

    OutputToGMSH("mfem_coarse_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
    std::cout << "  Coarse solution written to: mfem_coarse_solution.msh" << std::endl;

    // ========================================================================
    // Output fine-scale solution (MFEM reconstruction)
    // ========================================================================

    std::cout << "\nReconstructing fine-scale solution..." << std::endl;
    timer.reset();
    scalarLap.OutputMFEMFine(fullSolution.data(), dParams.numElePerDir[0], dParams.numElePerDir[1]);
    double reconstructTime = timer.seconds();
    std::cout << "  Fine-scale reconstruction time: " << reconstructTime * 1000.0 << " ms" << std::endl;
    std::cout << "  Fine solution written to: outputFine.txt" << std::endl;

    // ========================================================================
    // Performance summary
    // ========================================================================

    std::cout << "\n=== Performance Summary ===" << std::endl;
    std::cout << "MFEM Assembly time:       " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << "Factorization time:       " << factorTime * 1000.0 << " ms" << std::endl;
    std::cout << "Solve time:               " << solveTime * 1000.0 << " ms" << std::endl;
    std::cout << "Fine reconstruction time: " << reconstructTime * 1000.0 << " ms" << std::endl;
    std::cout << "Total time:               " << (assemblyTime + factorTime + solveTime + reconstructTime) * 1000.0 << " ms"
              << std::endl;
    std::cout << std::endl;

    // Compute some solution statistics
    double minVal = *std::min_element(solution.begin(), solution.end());
    double maxVal = *std::max_element(solution.begin(), solution.end());
    double avgVal = std::accumulate(solution.begin(), solution.end(), 0.0) / solution.size();

    std::cout << "Coarse solution statistics:" << std::endl;
    std::cout << "  Min value: " << minVal << std::endl;
    std::cout << "  Max value: " << maxVal << std::endl;
    std::cout << "  Avg value: " << avgVal << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
