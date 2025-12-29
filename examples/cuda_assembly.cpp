///
/// @file cuda_assembly.cpp
/// @brief Example demonstrating CUDA-accelerated FEM assembly
///
/// This example shows how to use the ScaledLaplacian class to assemble
/// the scaled Laplacian operator on CUDA devices using Kokkos.
///
/// Solves: -∇·(α∇u) = f on a 2D rectangular domain
/// with homogeneous Dirichlet boundary conditions
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
    std::cout << "=== CUDA FEM Assembly Example ===" << std::endl;
    std::cout << "Kokkos execution space: " << typeid(Kokkos::DefaultExecutionSpace).name() << std::endl;
    std::cout << "Precision: " << (sizeof(Scalar) == 4 ? "FP32 (float)" : "FP64 (double)")
              << " (" << sizeof(Scalar) << " bytes)" << std::endl;
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
    // Parse command line arguments
    // ========================================================================

    // Default values
    int nx = 8;
    int ny = 8;
    IMSI::ElementType elementType = IMSI::ElementType::Q1;
    bool writeOutput = false;

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
      } else if (arg == "-output" || arg == "--output") {
        writeOutput = true;
      } else if (arg == "-h" || arg == "--help") {
        std::cout << "\nUsage: " << argv[0] << " [options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -nx <n>    Number of elements in x direction (default: 8)" << std::endl;
        std::cout << "  -ny <n>    Number of elements in y direction (default: 8)" << std::endl;
        std::cout << "  -q1        Use Q1 (bilinear) elements (default)" << std::endl;
        std::cout << "  -mfem      Use MFEM_L (multiscale) elements" << std::endl;
        std::cout << "  -output    Write solution files (default: off)" << std::endl;
        std::cout << "  -h, --help Show this help message" << std::endl;
        std::cout << "\nExample: " << argv[0] << " -nx 16 -ny 16 -mfem -output" << std::endl;
        Kokkos::finalize();
        return 0;
      }
    }

    // Print coefficient type
#if IMSI_COEFF == 1
    std::cout << "\nCoefficients: CONSTANT (Laplace equation)" << std::endl;
#elif IMSI_COEFF == 2
    std::cout << "\nCoefficients: VARYING (oscillatory problem)" << std::endl;
#endif

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
    // Create sparse matrix structure on host
    // ========================================================================

    std::cout << "\nAllocating sparse matrix structure..." << std::endl;

    auto const numNodes = mesh.NumberVertices();
    auto const numNonZeros = meshConn.n2n.entries.extent(0);

    std::cout << "  Matrix size: " << numNodes << " x " << numNodes << std::endl;
    std::cout << "  Non-zeros:   " << numNonZeros << std::endl;

    // Host views
    Kokkos::View<Scalar*, Kokkos::HostSpace> rhs_h("rhs_host", numNodes);
    Kokkos::View<size_t*, Kokkos::HostSpace> matRowPtr_h("matRowPtr_host", numNodes + 1);
    Kokkos::View<int*, Kokkos::HostSpace> matColIdx_h("matColIdx_host", numNonZeros);
    Kokkos::View<Scalar*, Kokkos::HostSpace> matValues_h("matValues_host", numNonZeros);

    // Copy graph structure from n2n
    for (size_t i = 0; i <= numNodes; ++i) {
      matRowPtr_h(i) = meshConn.n2n.row_map(i);
    }
    for (size_t i = 0; i < numNonZeros; ++i) {
      matColIdx_h(i) = meshConn.n2n.entries(i);
      matValues_h(i) = 0.0;
    }

    // Device views
    Kokkos::View<Scalar*, Kokkos::Cuda> rhs_d("rhs_device", numNodes);
    Kokkos::View<size_t*, Kokkos::Cuda> matRowPtr_d("matRowPtr_device", numNodes + 1);
    Kokkos::View<int*, Kokkos::Cuda> matColIdx_d("matColIdx_device", numNonZeros);
    Kokkos::View<Scalar*, Kokkos::Cuda> matValues_d("matValues_device", numNonZeros);

    // Copy structure to device
    Kokkos::deep_copy(matRowPtr_d, matRowPtr_h);
    Kokkos::deep_copy(matColIdx_d, matColIdx_h);
    Kokkos::deep_copy(rhs_d, 0.0);
    Kokkos::deep_copy(matValues_d, 0.0);

    // ========================================================================
    // Assembly on CUDA
    // ========================================================================

    std::cout << "\n=== Starting CUDA Assembly ===" << std::endl;

    // Create coefficient functor instances
    AxFunctor ax_func;
    AyFunctor ay_func;
    FFunctor f_func;

    // Instantiate ScaledLaplacian with ExecutionSpace and functor types
    auto scalarLap = ScaledLaplacian<Kokkos::Cuda, AxFunctor, AyFunctor, FFunctor>(
        meshConn, RuleType::Gauss, 2, ax_func, ay_func, f_func);

    Kokkos::Timer timer;
    if (elementType == IMSI::ElementType::MFEM_L) {
      scalarLap.GetLinearSystemMFEM(rhs_d, matRowPtr_d, matColIdx_d, matValues_d);
    } else {
      scalarLap.GetLinearSystem(rhs_d, matRowPtr_d, matColIdx_d, matValues_d);
    }
    Kokkos::fence();
    double assemblyTime = timer.seconds();

    std::cout << "=== Assembly Complete ===" << std::endl;
    std::cout << "Assembly time: " << assemblyTime * 1000.0 << " ms" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // Copy results back to host
    // ========================================================================

    timer.reset();
    Kokkos::deep_copy(rhs_h, rhs_d);
    Kokkos::deep_copy(matValues_h, matValues_d);
    Kokkos::fence();
    double copyTime = timer.seconds();
    std::cout << "Copying (Kc, fc) to host: " << copyTime * 1000.0 << " ms" << std::endl;

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

    std::vector<Scalar> fullSolution(numNodes, 0.0);
    for (size_t i = 0; i < numFreeDofs; ++i) {
      fullSolution[freeToGlobal[i]] = solution[i];
    }

    double reconstructTime = 0.0;
    if (writeOutput) {
      std::cout << "\nWriting solution to file..." << std::endl;

      if (elementType == IMSI::ElementType::MFEM_L) {
        OutputToGMSH("cuda_coarse_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
        std::cout << "  Coarse solution written to: cuda_coarse_solution.msh" << std::endl;

        // Reconstruct fine-scale solution
        std::cout << "\nReconstructing fine-scale solution..." << std::endl;
        timer.reset();
        scalarLap.OutputMFEMFine(fullSolution.data(), nx, ny);
        reconstructTime = timer.seconds();
        std::cout << "  Fine-scale reconstruction time: " << reconstructTime * 1000.0 << " ms" << std::endl;
        std::cout << "  Fine solution written to: outputFine.txt" << std::endl;
      } else {
        OutputToGMSH("cuda_solution.msh", mesh, fullSolution.data(), int(fullSolution.size()));
        std::cout << "  Solution written to: cuda_solution.msh" << std::endl;
      }
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
