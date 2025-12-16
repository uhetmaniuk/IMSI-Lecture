///
/// @file openmp_assembly.cpp
/// @brief Example demonstrating OpenMP-accelerated FEM assembly using ScaledLaplacian
///
/// This example shows how to use the ScaledLaplacian class with Kokkos::OpenMP
/// execution space to assemble the scaled Laplacian operator using OpenMP threads.
///
/// Solves: -∇·(α∇u) = f on a 2D rectangular domain
/// with homogeneous Dirichlet boundary conditions
///
/// NOTE: ScaledLaplacian works with any Kokkos execution space including
///       OpenMP, CUDA, HIP, etc. due to its functor-based design pattern.
///

#include <Kokkos_Core.hpp>
#include <iostream>
#include <vector>

// Defines Scalar, AxFunctor, AyFunctor, FFunctor
#include "config_assembly.h"

#include "../src/Element.h"
#include "../src/Mesh.h"
#include "../src/MeshUtils.h"
#include "../src/ScaledLaplacian.h"
#include "../src/SymmetricSparse.hpp"
#include "../src/Utils.h"

using namespace IMSI;

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::cout << "=== OpenMP FEM Assembly Example (using ScaledLaplacian) ===" << std::endl;
    std::cout << "Kokkos execution space: " << typeid(Kokkos::OpenMP).name() << std::endl;
    std::cout << "Precision: " << (sizeof(Scalar) == 4 ? "FP32 (float)" : "FP64 (double)")
              << " (" << sizeof(Scalar) << " bytes)" << std::endl;
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
    // Parse command line arguments
    // ========================================================================

    // Default values
    int nx = 256;
    int ny = 256;
    IMSI::ElementType elementType = IMSI::ElementType::Q1;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "-nx" && i + 1 < argc) {
        nx = std::atoi(argv[++i]);
      } else if (arg == "-ny" && i + 1 < argc) {
        ny = std::atoi(argv[++i]);
      } else if (arg == "-q1") {
        elementType = IMSI::ElementType::Q1;
      } else if (arg == "-mfem") {
        elementType = IMSI::ElementType::MFEM_L;
      } else if (arg == "-h" || arg == "--help") {
        std::cout << "\nUsage: " << argv[0] << " [options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -nx <n>    Number of elements in x direction (default: 256)" << std::endl;
        std::cout << "  -ny <n>    Number of elements in y direction (default: 256)" << std::endl;
        std::cout << "  -q1        Use Q1 (bilinear) elements (default)" << std::endl;
        std::cout << "  -mfem      Use MFEM_L (multiscale) elements" << std::endl;
        std::cout << "  -h, --help Show this help message" << std::endl;
        std::cout << "\nExample: " << argv[0] << " -nx 16 -ny 16 -q1" << std::endl;
        Kokkos::finalize();
        return 0;
      }
    }

    // ========================================================================
    // Problem setup
    // ========================================================================

    // ========================================================================
    // Mesh generation
    // ========================================================================

    IMSI::DomainParams dParams;
    dParams.numElePerDir[0] = nx;
    dParams.numElePerDir[1] = ny;
    dParams.omega           = IMSI::DomainType::Rectangle;
    dParams.cellType        = elementType;

    std::string elementName = (elementType == IMSI::ElementType::Q1) ? "Q1" : "MFEM_L";
    std::cout << "\nGenerating mesh: " << dParams.numElePerDir[0] << " x " << dParams.numElePerDir[1]
              << " " << elementName << " elements" << std::endl;

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
    // Create sparse matrix structure
    // ========================================================================

    std::cout << "\nAllocating sparse matrix structure..." << std::endl;

    auto const numNodes = mesh.NumberVertices();
    auto const numNonZeros = meshConn.n2n.entries.extent(0);

    std::cout << "  Matrix size: " << numNodes << " x " << numNodes << std::endl;
    std::cout << "  Non-zeros:   " << numNonZeros << std::endl;

    // OpenMP execution space views (host-accessible)
    Kokkos::View<Scalar*, Kokkos::OpenMP> rhs("rhs", numNodes);
    Kokkos::View<size_t*, Kokkos::OpenMP> matRowPtr("matRowPtr", numNodes + 1);
    Kokkos::View<int*, Kokkos::OpenMP> matColIdx("matColIdx", numNonZeros);
    Kokkos::View<Scalar*, Kokkos::OpenMP> matValues("matValues", numNonZeros);

    // Copy graph structure from n2n (on host)
    for (size_t i = 0; i <= numNodes; ++i) {
      matRowPtr(i) = meshConn.n2n.row_map(i);
    }
    for (size_t i = 0; i < numNonZeros; ++i) {
      matColIdx(i) = meshConn.n2n.entries(i);
      matValues(i) = Scalar(0);
    }

    // Initialize to zero
    Kokkos::deep_copy(rhs, Scalar(0));

    // ========================================================================
    // Assembly with OpenMP
    // ========================================================================

    std::cout << "\n=== Starting OpenMP Assembly ===" << std::endl;

    // Create coefficient functor instances
    AxFunctor ax_func;
    AyFunctor ay_func;
    FFunctor f_func;

    // Instantiate ScaledLaplacian with OpenMP ExecutionSpace and functor types
    auto scalarLap = ScaledLaplacian<Kokkos::OpenMP, AxFunctor, AyFunctor, FFunctor>(
        meshConn, RuleType::Gauss, 2, ax_func, ay_func, f_func);

    Kokkos::Timer timer;
    if (elementType == IMSI::ElementType::MFEM_L) {
      scalarLap.GetLinearSystemMFEM(rhs, matRowPtr, matColIdx, matValues);
    } else {
      scalarLap.GetLinearSystem(rhs, matRowPtr, matColIdx, matValues);
    }
    Kokkos::fence();
    double assemblyTime = timer.seconds();

    std::cout << "=== Assembly Complete ===" << std::endl;
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
    std::cout << "  Total DOFs:    " << numNodes << std::endl;
    std::cout << "  Boundary DOFs: " << numNodes - numFreeDofs << std::endl;
    std::cout << "  Free DOFs:     " << numFreeDofs << std::endl;

    // Build reduced system
    std::vector<int> newRowPtr(numFreeDofs + 1, 0);
    std::vector<int> newColIdx;
    std::vector<Scalar> newValues;
    std::vector<Scalar> newRhs(numFreeDofs);

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

    SymmetricSparse<Scalar> solver(
        numFreeDofs, newColIdx.size(), newRowPtr.data(), newColIdx.data(), newValues.data(), true);

    timer.reset();
    solver.factor();
    double factorTime = timer.seconds();
    std::cout << "  Factorization time: " << factorTime * 1000.0 << " ms" << std::endl;

    std::vector<Scalar> solution(numFreeDofs);
    timer.reset();
    solver.Solve(1, newRhs.data(), solution.data());
    double solveTime = timer.seconds();
    std::cout << "  Solve time:         " << solveTime * 1000.0 << " ms" << std::endl;

    // ========================================================================
    // Output results
    // ========================================================================

    std::cout << "\nWriting solution to file..." << std::endl;

    std::vector<Scalar> fullSolution(numNodes, 0.0);
    for (size_t i = 0; i < numFreeDofs; ++i) {
      fullSolution[freeToGlobal[i]] = solution[i];
    }

    double reconstructTime = 0.0;
    if (elementType == IMSI::ElementType::MFEM_L) {
      OutputToGMSH("openmp_coarse_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
      std::cout << "  Coarse solution written to: openmp_coarse_solution.msh" << std::endl;

      // Reconstruct fine-scale solution
      std::cout << "\nReconstructing fine-scale solution..." << std::endl;
      timer.reset();
      scalarLap.OutputMFEMFine(fullSolution.data(), nx, ny);
      reconstructTime = timer.seconds();
      std::cout << "  Fine-scale reconstruction time: " << reconstructTime * 1000.0 << " ms" << std::endl;
      std::cout << "  Fine solution written to: outputFine.txt" << std::endl;
    } else {
      OutputToGMSH("openmp_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
      std::cout << "  Solution written to: openmp_solution.msh" << std::endl;
    }

    // ========================================================================
    // Performance summary
    // ========================================================================

    std::cout << "\n=== Performance Summary ===" << std::endl;
    std::cout << "Assembly time:       " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << "Factorization time:  " << factorTime * 1000.0 << " ms" << std::endl;
    std::cout << "Solve time:          " << solveTime * 1000.0 << " ms" << std::endl;
    if (elementType == IMSI::ElementType::MFEM_L) {
      std::cout << "Reconstruction time: " << reconstructTime * 1000.0 << " ms" << std::endl;
      std::cout << "Total time:          " << (assemblyTime + factorTime + solveTime + reconstructTime) * 1000.0 << " ms"
                << std::endl;
    } else {
      std::cout << "Total time:          " << (assemblyTime + factorTime + solveTime) * 1000.0 << " ms"
                << std::endl;
    }
    std::cout << std::endl;

    // Compute some solution statistics
    Scalar minVal = *std::min_element(solution.begin(), solution.end());
    Scalar maxVal = *std::max_element(solution.begin(), solution.end());
    Scalar avgVal = std::accumulate(solution.begin(), solution.end(), Scalar(0.0)) / solution.size();

    std::cout << "Solution statistics:" << std::endl;
    std::cout << "  Min value: " << minVal << std::endl;
    std::cout << "  Max value: " << maxVal << std::endl;
    std::cout << "  Avg value: " << avgVal << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
