///
/// @file cuda_assembly.cpp
/// @brief Example demonstrating CUDA-accelerated FEM assembly
///
/// This example shows how to use the ScaledLaplacianCuda class to assemble
/// the scaled Laplacian operator on CUDA devices using Kokkos.
///
/// Solves: -∇·(α∇u) = f on a 2D rectangular domain
/// with homogeneous Dirichlet boundary conditions
///

#include <Kokkos_Core.hpp>
#include <cmath>
#include <functional>
#include <iostream>

#include "../src/Element.h"
#include "../src/Mesh.h"
#include "../src/MeshUtils.h"
#include "../src/ScaledLaplacianCuda.h"
#include "../src/SymmetricSparse.hpp"
#include "../src/Utils.h"

using namespace IMSI;

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::cout << "=== CUDA FEM Assembly Example ===" << std::endl;
    std::cout << "Kokkos execution space: " << typeid(Kokkos::DefaultExecutionSpace).name() << std::endl;
    std::cout << std::endl;

    // Check if CUDA is available
#ifdef KOKKOS_ENABLE_CUDA
    std::cout << "CUDA is enabled" << std::endl;
    std::cout << "CUDA device count: " << Kokkos::Cuda().concurrency() << " threads" << std::endl;
#else
    std::cerr << "ERROR: This example requires Kokkos with CUDA support" << std::endl;
    std::cerr << "       Rebuild with -DKokkos_ENABLE_CUDA=ON" << std::endl;
    Kokkos::finalize();
    return 1;
#endif

    // ========================================================================
    // Problem setup
    // ========================================================================

    // Define material coefficients and forcing term
    auto alpha_x = [](double x, double y, double z) -> double { return 1.0; };
    auto beta_y  = [](double x, double y, double z) -> double { return 1.0; };
    auto f_rhs   = [](double x, double y, double z) -> double {
      // Manufactured solution: u = sin(pi*x) * sin(pi*y)
      // => -Laplacian(u) = 2*pi^2 * sin(pi*x) * sin(pi*y)
      constexpr double pi = 3.14159265358979323846;
      return 2.0 * pi * pi * std::sin(pi * x) * std::sin(pi * y);
    };

    // ========================================================================
    // Mesh generation
    // ========================================================================

    IMSI::DomainParams dParams;
    dParams.numElePerDir[0] = 64; // 2048
    dParams.numElePerDir[1] = 64; // 2048
    dParams.omega           = IMSI::DomainType::Rectangle;
    dParams.cellType        = IMSI::ElementType::MFEM_L;

    std::cout << "Generating mesh: " << dParams.numElePerDir[0] << " x " << dParams.numElePerDir[1]
              << " Q1 elements" << std::endl;

    auto mesh = IMSI::GenerateMesh(dParams);
    std::cout << "  Number of elements: " << mesh.NumberCells() << std::endl;
    std::cout << "  Number of nodes:    " << mesh.NumberVertices() << std::endl;

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
    // Create sparse matrix structure on host
    // ========================================================================

    std::cout << "\nAllocating sparse matrix structure..." << std::endl;

    auto const numNodes = mesh.NumberVertices();
    auto const numNonZeros = meshConn.n2n.entries.extent(0);

    std::cout << "  Matrix size: " << numNodes << " x " << numNodes << std::endl;
    std::cout << "  Non-zeros:   " << numNonZeros << std::endl;

    // Host views
    Kokkos::View<double*, Kokkos::HostSpace> rhs_h("rhs_host", numNodes);
    Kokkos::View<size_t*, Kokkos::HostSpace> matRowPtr_h("matRowPtr_host", numNodes + 1);
    Kokkos::View<int*, Kokkos::HostSpace> matColIdx_h("matColIdx_host", numNonZeros);
    Kokkos::View<double*, Kokkos::HostSpace> matValues_h("matValues_host", numNonZeros);

    // Copy graph structure from n2n
    for (size_t i = 0; i <= numNodes; ++i) {
      matRowPtr_h(i) = meshConn.n2n.row_map(i);
    }
    for (size_t i = 0; i < numNonZeros; ++i) {
      matColIdx_h(i) = meshConn.n2n.entries(i);
      matValues_h(i) = 0.0;
    }

    // Device views
    Kokkos::View<double*, Kokkos::Cuda> rhs_d("rhs_device", numNodes);
    Kokkos::View<size_t*, Kokkos::Cuda> matRowPtr_d("matRowPtr_device", numNodes + 1);
    Kokkos::View<int*, Kokkos::Cuda> matColIdx_d("matColIdx_device", numNonZeros);
    Kokkos::View<double*, Kokkos::Cuda> matValues_d("matValues_device", numNonZeros);

    // Copy structure to device
    Kokkos::deep_copy(matRowPtr_d, matRowPtr_h);
    Kokkos::deep_copy(matColIdx_d, matColIdx_h);
    Kokkos::deep_copy(rhs_d, 0.0);
    Kokkos::deep_copy(matValues_d, 0.0);

    // ========================================================================
    // Assembly on CUDA
    // ========================================================================

    std::cout << "\n=== Starting CUDA Assembly ===" << std::endl;

    ScaledLaplacianCuda scalarLap(meshConn, alpha_x, beta_y, f_rhs, RuleType::Gauss, 2);

    Kokkos::Timer timer;
    scalarLap.GetLinearSystem(rhs_d, matRowPtr_d, matColIdx_d, matValues_d);
    Kokkos::fence();
    double assemblyTime = timer.seconds();

    std::cout << "=== Assembly Complete ===" << std::endl;
    std::cout << "Assembly time: " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // Copy results back to host
    // ========================================================================

    std::cout << "Copying results to host..." << std::endl;
    Kokkos::deep_copy(rhs_h, rhs_d);
    Kokkos::deep_copy(matValues_h, matValues_d);

    // ========================================================================
    // Apply boundary conditions
    // ========================================================================

    std::cout << "Applying boundary conditions..." << std::endl;

    auto const&      bdyNodes = myMesh.GetBoundaryNodes();
    std::vector<int> globalToFree(myMesh.NumberVertices());
    std::vector<int> freeToGlobal(size(globalToFree) - size(bdyNodes));
    IMSI::MapDegreesOfFreedom(bdyNodes, globalToFree, freeToGlobal);

    std::cout << "  Total DOFs:    " << numNodes << std::endl;
    std::cout << "  Boundary DOFs: " << numNodes - numFreeDofs << std::endl;
    std::cout << "  Free DOFs:     " << numFreeDofs << std::endl;

    // Build reduced system
    std::vector<int> newRowPtr(numFreeDofs + 1, 0);
    std::vector<int> newColIdx;
    std::vector<double> newValues;
    std::vector<double> newRhs(numFreeDofs);

    newRowPtr[0] = 0;
    for (size_t i = 0; i < numFreeDofs; ++i) {
      auto gDof = freeToGlobal[i];
      newRhs[i] = rhs_h(gDof);
      for (size_t k = matRowPtr_h(gDof); k < matRowPtr_h(gDof + 1); ++k) {
        auto gCol = matColIdx_h(k);
        if (globalToFree[gCol] != -1) {
          newColIdx.push_back(globalToFree[gCol]);
          newValues.push_back(matValues_h(k));
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
    // Output results
    // ========================================================================

    std::cout << "\nWriting solution to file..." << std::endl;

    std::vector<double> fullSolution(numNodes, 0.0);
    for (size_t i = 0; i < numFreeDofs; ++i) {
      fullSolution[freeToGlobal[i]] = solution[i];
    }

    OutputToGMSH("cuda_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
    std::cout << "  Solution written to: cuda_solution.msh" << std::endl;

    // ========================================================================
    // Performance summary
    // ========================================================================

    std::cout << "\n=== Performance Summary ===" << std::endl;
    std::cout << "Assembly time:       " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << "Factorization time:  " << factorTime * 1000.0 << " ms" << std::endl;
    std::cout << "Solve time:          " << solveTime * 1000.0 << " ms" << std::endl;
    std::cout << "Total time:          " << (assemblyTime + factorTime + solveTime) * 1000.0 << " ms"
              << std::endl;
    std::cout << std::endl;

    // Compute some solution statistics
    double minVal = *std::min_element(solution.begin(), solution.end());
    double maxVal = *std::max_element(solution.begin(), solution.end());
    double avgVal = std::accumulate(solution.begin(), solution.end(), 0.0) / solution.size();

    std::cout << "Solution statistics:" << std::endl;
    std::cout << "  Min value: " << minVal << std::endl;
    std::cout << "  Max value: " << maxVal << std::endl;
    std::cout << "  Avg value: " << avgVal << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
