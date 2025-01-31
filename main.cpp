#include <chrono>
#include <iostream>
#include <vector>

#include "Kokkos_Core.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"
#include "src/FunctionExamples.h"
#include "src/Mesh.h"
#include "src/MeshUtils.h"
#include "src/ScaledLaplacian.h"
#include "src/Utils.h"

using accelerator_space = Kokkos::DefaultExecutionSpace;
using accelerator_type  = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;

using execution_space     = Kokkos::DefaultHostExecutionSpace;
using host_execution_type = typename Tacho::UseThisDevice<execution_space>::type;

int
main(int argc, char* argv[])
{
  //
  // Flag --kokkos-num-threads=INT
  //
  Kokkos::initialize(argc, argv);
  {
    //
    std::cout << " ## THREADS " << Kokkos::num_threads() << "\n";
    //
    IMSI::DomainParams dParams;
    dParams.numElePerDir[0] = 1024;
    dParams.numElePerDir[1] = 1024;
    //
    std::cout << " Grid = " << dParams.numElePerDir[0] << " x " << dParams.numElePerDir[1] << "\n";
    //
    auto                          start  = std::chrono::high_resolution_clock::now();
    auto                          myMesh = IMSI::GenerateMesh(dParams);
    auto                          end    = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt     = end - start;
    std::cout << " --- Generate Mesh = " << dt.count() << "\n";
    std::cout << " Vertices = " << myMesh.NumberVertices() << "\n";
    //
    start         = std::chrono::high_resolution_clock::now();
    auto meshData = IMSI::GetMeshConnectivity(myMesh);
    end           = std::chrono::high_resolution_clock::now();
    dt            = end - start;
    std::cout << " --- Get Mesh Connectivities = " << dt.count() << "\n";
    //
    auto const                             numDofs = meshData.mesh.NumberVertices();
    Kokkos::View<size_t*, execution_space> matRowPtr("Row Pointer Matrix", numDofs + 1);
    Kokkos::deep_copy(matRowPtr, meshData.n2n.row_map);
    Kokkos::View<int*, execution_space> matColIdx("Column Index Matrix", meshData.n2n.entries.size());
    Kokkos::deep_copy(matColIdx, meshData.n2n.entries);
    Kokkos::View<double*, execution_space> matValues("Values Matrix", matColIdx.size());
    Kokkos::View<double*, execution_space> rhsValues("Values RHS", numDofs);
    //
    const IMSI::ParabolicPb problem;
    start = std::chrono::high_resolution_clock::now();
    IMSI::ScaledLaplacian dataAssembly(meshData, problem.ax, problem.ay, problem.f, IMSI::RuleType::Gauss, 2);
    dataAssembly.GetLinearSystem(rhsValues, matRowPtr, matColIdx, matValues);
    end = std::chrono::high_resolution_clock::now();
    dt  = end - start;
    std::cout << " --- GetLinearSystem = " << dt.count() << "\n";
    //
    //--- Get map from global numbering to free numbering (and vice versa)
    //
    auto const&      bdyNodes = myMesh.GetBoundaryNodes();
    std::vector<int> globalToFree(myMesh.NumberVertices());
    std::vector<int> freeToGlobal(size(globalToFree) - size(bdyNodes));
    IMSI::MapDegreesOfFreedom(bdyNodes, globalToFree, freeToGlobal);
    //
    //--- Solve the linear system
    //
    Kokkos::View<double*, execution_space> u("Value for approximation", myMesh.NumberVertices());
    {
      Kokkos::View<size_t*, execution_space> newRowPtr("Row Pointer for free degrees", size(freeToGlobal) + 1);
      int                                    newNNZ = 0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, size(freeToGlobal)),
          KOKKOS_LAMBDA(int i, int& count) {
            auto   gDof   = freeToGlobal[i];
            size_t iCount = 0;
            for (auto k = matRowPtr[gDof]; k < matRowPtr[gDof + 1]; ++k) {
              iCount += (globalToFree[matColIdx[k]] != -1);
            }
            newRowPtr[i + 1] = iCount;
            count += iCount;
          },
          Kokkos::Sum<int>(newNNZ));
      //
      auto h_rowPtr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), newRowPtr);
      for (int i = 0; i < h_rowPtr.size() - 1; ++i) { h_rowPtr(i + 1) += h_rowPtr(i); }
      Kokkos::deep_copy(newRowPtr, h_rowPtr);
      //
      Kokkos::View<int*, execution_space>    newColIdx("Col Idx for free degrees", newNNZ);
      Kokkos::View<double*, execution_space> newValues("Values for free degrees", newNNZ);
      Kokkos::View<double*, execution_space> newRHS("RHS for free degrees", size(freeToGlobal));
      //
      Kokkos::parallel_for(
          "Fill free Matrix", Kokkos::RangePolicy<execution_space>(0, size(freeToGlobal)), KOKKOS_LAMBDA(int iFree) {
            auto   gdof = freeToGlobal[iFree];
            size_t pos  = newRowPtr[iFree];
            for (auto k = matRowPtr[gdof]; k < matRowPtr[gdof + 1]; ++k) {
              auto const gCol = matColIdx[k];
              if (globalToFree[gCol] != -1) {
                newColIdx[pos] = globalToFree[gCol];
                newValues[pos] = matValues[k];
                pos += 1;
              }
            }
            newRHS[iFree] = rhsValues[gdof];
          });
      //
      Tacho::CrsMatrixBase<double, host_execution_type> h_A;
      auto const                                        n = size(freeToGlobal);
      h_A.setExternalMatrix(n, n, newNNZ, newRowPtr, newColIdx, newValues);
      //
      Tacho::CrsMatrixBase<double, accelerator_type> A;
      A.createMirror(h_A);
      Tacho::Driver<double, accelerator_type> solver;
      solver.analyze(A.NumRows(), A.RowPtr(), A.Cols());
      //
      Kokkos::Timer timer;
      solver.initialize();
      double initi_time = timer.seconds();
      timer.reset();
      solver.factorize(A.Values());
      double facto_time = timer.seconds();
      //
      // Tacho::Driver::solve expects a "matrix" for parameters.
      // So even if one linear system with a single RHS needs to be solved,
      // the current interface requires to use a "matrix" with 1 column.
      //
      Kokkos::View<double**, Kokkos::LayoutLeft, typename accelerator_type::execution_space> b("rhs", n, 1);
      Kokkos::View<double**, Kokkos::LayoutLeft, typename accelerator_type::execution_space> x("solution", n, 1);
      Kokkos::View<double**, Kokkos::LayoutLeft, typename accelerator_type::execution_space> wt("workspace", n, 1);
      //
      // Fill the right hand side with random values
      auto b0 = Kokkos::subview(b, Kokkos::ALL, 0);
      Kokkos::deep_copy(b0, newRHS);
      //
      double solve_time = 0.0;
      timer.reset();
      solver.solve(x, b, wt);
      solve_time += timer.seconds();
      {
        const double res = solver.computeRelativeResidual(A.Values(), x, b);
        std::cout << "TachoSolver: residual = " << res << "\n";
      }
      std::cout << std::endl;
      std::cout << " Initialize Time " << initi_time << std::endl;
      std::cout << " Facto Time " << facto_time << std::endl;
      std::cout << " Solve Time " << solve_time << std::endl;
      std::cout << std::endl;
      solver.release();
      //
      auto ufree = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(x, Kokkos::ALL, 0));
      Kokkos::parallel_for(
          "Fill u", Kokkos::RangePolicy<execution_space>(0, size(freeToGlobal)), KOKKOS_LAMBDA(int i) {
            auto dof = freeToGlobal[i];
            u[dof]   = ufree[i];
          });
    }
    //
    //--- Compute error norms
    //
    //
    //--- Output the approximation to GMSH
    //
    IMSI::OutputToGMSH("approximate.msh", myMesh, u.data(), myMesh.NumberVertices());
    //
  }
  Kokkos::finalize();
  return 0;
}
