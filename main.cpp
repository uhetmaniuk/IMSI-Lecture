#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#ifdef IMSI_USE_TACHO
#include "TPL/tacho/src/Tacho_CrsMatrixBase.hpp"
#include "TPL/tacho/src/Tacho_Driver.hpp"
#endif

#include "main_config.h"
#include "src/FunctionExamples.h"
#include "src/Mesh.h"
#include "src/MeshUtils.h"
#include "src/ScaledLaplacian.h"
#include "src/SymmetricSparse.hpp"
#include "src/Utils.h"

int
main(int argc, char* argv[])
{
  //
  // Flag --kokkos-num-threads=INT
  //
  bool useSIMD     = false;
  bool useColoring = true;
  //
  for (int ii = 0; ii < argc; ++ii) {
    std::string arg = argv[ii];
    if (arg.compare("--simd") == 0) { useSIMD = true; }
    if (arg.compare("--nocolor") == 0) { useColoring = false; }
  }
  Kokkos::initialize(argc, argv);
  {
    //
    std::cout << " ## THREADS " << Kokkos::num_threads() << "\n";
    //
    {
      /*
      /// HACK
      typedef Kokkos::TeamPolicy<accelerator_space>::member_type member_type;
      // Launch a kernel
      auto uhStart = std::chrono::high_resolution_clock::now();
      Kokkos::parallel_for("Test RP", Kokkos::TeamPolicy<accelerator_space>(4, Kokkos::AUTO, 16),
      KOKKOS_LAMBDA(member_type team_member) {
        /// league_rank -> Colored element ID
          /// team -> only 1 team
          ///
        printf(" lsize %d lrank %d tsize %d trank %d \n", team_member.league_size(), team_member.league_rank(),
      team_member.team_size(), team_member.team_rank()); std::this_thread::sleep_for(std::chrono::seconds(3));
        Kokkos::parallel_for (Kokkos::TeamThreadRange(team_member, 7),
            KOKKOS_LAMBDA (const int& i)
        {
          printf(" TeamThreadRange -> %d \n", i);
          std::this_thread::sleep_for(std::chrono::seconds(1));
        });
        Kokkos::parallel_for (Kokkos::TeamVectorRange(team_member, 11),
            KOKKOS_LAMBDA (const int& i)
        {
          printf(" TeamVectorRange -> %d \n", i);
          std::this_thread::sleep_for(std::chrono::seconds(2));
        });
      });
      auto uhStop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> uhDT = uhStop - uhStart;
      printf(" dt %e \n", uhDT.count());
      */
    }
    //
    ///
    IMSI::DomainParams dParams;
    dParams.numElePerDir[0] = 8; // 2048
    dParams.numElePerDir[1] = 8; // 2048
    dParams.omega           = IMSI::DomainType::Rectangle;
    dParams.cellType        = IMSI::ElementType::Q1;
    //
    std::cout << " Grid = " << dParams.numElePerDir[0] << " x " << dParams.numElePerDir[1] << "\n";
    if (dParams.cellType == IMSI::ElementType::MFEM_L) { std::cout << " MFEM_L (-Q1) discretization \n"; }
    if (dParams.cellType == IMSI::ElementType::Q1) { std::cout << " Q1 discretization \n"; }
    if (dParams.cellType == IMSI::ElementType::Q2) { std::cout << " Q2 discretization \n"; }
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
    auto const                                  numDofs = meshData.mesh.NumberVertices();
    Kokkos::View<size_t*, host_execution_space> matRowPtr("Row Pointer Matrix", numDofs + 1);
    Kokkos::deep_copy(matRowPtr, meshData.n2n.row_map);
    Kokkos::View<int*, host_execution_space> matColIdx("Column Index Matrix", meshData.n2n.entries.size());
    Kokkos::deep_copy(matColIdx, meshData.n2n.entries);
    Kokkos::View<double*, host_execution_space> matValues("Values Matrix", matColIdx.size());
    Kokkos::View<double*, host_execution_space> rhsValues("Values RHS", numDofs);
    //
    const IMSI::ParabolicPb2 problem;
    start = std::chrono::high_resolution_clock::now();
    IMSI::ScaledLaplacian dataAssembly(meshData, problem.ax, problem.ay, problem.f, IMSI::RuleType::Gauss, 4);
    if (useColoring) {
      if ((useSIMD) && (dParams.cellType != IMSI::ElementType::MFEM_L)) {
        dataAssembly.GetLinearSystem_v<host_execution_space, true>(rhsValues, matRowPtr, matColIdx, matValues);
      } else {
        dataAssembly.GetLinearSystem_v<host_execution_space, false>(rhsValues, matRowPtr, matColIdx, matValues);
      }
    } else {
      if ((useSIMD) && (dParams.cellType != IMSI::ElementType::MFEM_L)) {
        dataAssembly.GetLinearSystem_v<host_execution_space, true, false>(rhsValues, matRowPtr, matColIdx, matValues);
      } else {
        dataAssembly.GetLinearSystem_v<host_execution_space, false, false>(rhsValues, matRowPtr, matColIdx, matValues);
      }
    }
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
    Kokkos::View<double*, host_execution_space> u("Value for approximation", myMesh.NumberVertices());
    {
      Kokkos::View<int*, host_execution_space> newRowPtr("Row Pointer for free degrees", size(freeToGlobal) + 1);
      int                                      newNNZ = 0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<host_execution_space>(0, size(freeToGlobal)),
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
      for (int i = 0; i < newRowPtr.size() - 1; ++i) { newRowPtr(i + 1) += newRowPtr(i); }
      //
      Kokkos::View<int*, host_execution_space>    newColIdx("Col Idx for free degrees", newNNZ);
      Kokkos::View<double*, host_execution_space> newValues("Values for free degrees", newNNZ);
      Kokkos::View<double*, host_execution_space> newRHS("RHS for free degrees", size(freeToGlobal));
      //
      auto const n = size(freeToGlobal);
      Kokkos::parallel_for(
          "Fill free Matrix", Kokkos::RangePolicy<host_execution_space>(0, n), KOKKOS_LAMBDA(int iFree) {
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
      Kokkos::View<double*, host_execution_space> uFree("Values for free degrees", n);
      SymmetricSparse<double> Ktmp(n, newNNZ, newRowPtr.data(), newColIdx.data(), newValues.data(), true);
      Kokkos::Timer           timer;
      Ktmp.factor();
      double facto_time = timer.seconds();
      timer.reset();
      Ktmp.Solve(1, newRHS.data(), uFree.data());
      double solve_time = timer.seconds();
      //
      /*
      Tacho::CrsMatrixBase<double, host_execution_type> h_A;
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
      Kokkos::View<double**, Kokkos::LayoutLeft, accelerator_space> b("rhs", n, 1);
      Kokkos::View<double**, Kokkos::LayoutLeft, accelerator_space> x("solution", n, 1);
      Kokkos::View<double**, Kokkos::LayoutLeft, accelerator_space> wt("workspace", n, 1);
      //
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
      */
      std::cout << std::endl;
      std::cout << " Facto Time " << facto_time << std::endl;
      std::cout << " Solve Time " << solve_time << std::endl;
      std::cout << std::endl;
      //
      //      auto ufree = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(x, Kokkos::ALL, 0));
      //      for (int ii = 0; ii < std::min<int>(n, 128); ++ii) {
      //        printf(" u %e z %e diff %e \n", ufree[ii], z[ii], ufree[ii] - z[ii]);
      //      }
      Kokkos::parallel_for(
          "Fill u", Kokkos::RangePolicy<host_execution_space>(0, n), KOKKOS_LAMBDA(int i) {
            auto dof = freeToGlobal[i];
            u[dof]   = uFree[i];
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
    if (dParams.cellType == IMSI::ElementType::MFEM_L) {
      dataAssembly.OutputMFEMFine(u.data(), dParams.numElePerDir[0], dParams.numElePerDir[1]);
    }
  }
  Kokkos::finalize();
  return 0;
}
