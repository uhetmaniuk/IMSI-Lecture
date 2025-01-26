#include <chrono>
#include <iostream>
#include <vector>

#include "src/FunctionExamples.h"
#include "src/Mesh.h"
#include "src/MeshUtils.h"
#include "src/ScaledLaplacian.h"
#include "src/Utils.h"

#include "Kokkos_Core.hpp"

int main(int argc, char *argv[]) {
    //
    // Flag --kokkos-num-threads=INT
    //
    Kokkos::initialize(argc, argv);
    {
        using execution_space = typename Kokkos::DefaultHostExecutionSpace;
        //
        std::cout << " ## THREADS " << Kokkos::num_threads() << "\n";
        //
        IMSI::DomainParams dParams;
        dParams.numElePerDir[0] = 6;
        dParams.numElePerDir[1] = 7;
        //
        auto start = std::chrono::high_resolution_clock::now();
        auto myMesh = IMSI::GenerateMesh(dParams);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = end - start;
        std::cout << " --- Generate Mesh = " << dt.count() << "\n";
        //
        start = std::chrono::high_resolution_clock::now();
        auto meshData = IMSI::GetMeshConnectivity(myMesh);
        end = std::chrono::high_resolution_clock::now();
        dt = end - start;
        std::cout << " --- Get Mesh Connectivities = " << dt.count() << "\n";
        //
        auto const numDofs = meshData.mesh.NumberVertices();
        Kokkos::View<size_t *, execution_space> matRowPtr("Row Pointer Matrix", numDofs + 1);
        Kokkos::deep_copy(matRowPtr, meshData.n2n.row_map);
        Kokkos::View<int *, execution_space> matColIdx("Column Index Matrix", meshData.n2n.entries.size());
        Kokkos::deep_copy(matColIdx, meshData.n2n.entries);
        Kokkos::View<double *, execution_space> matValues("Values Matrix", matColIdx.size());
        Kokkos::View<double *, execution_space> rhsValues("Values RHS", numDofs);
        //
        const IMSI::ParabolicPb problem;
        start = std::chrono::high_resolution_clock::now();
        IMSI::ScaledLaplacian dataAssembly(meshData, problem.ax, problem.ay, problem.f, IMSI::RuleType::Gauss, 2);
        dataAssembly.GetLinearSystem(rhsValues, matRowPtr, matColIdx, matValues);
        end = std::chrono::high_resolution_clock::now();
        dt = end - start;
        std::cout << " --- GetLinearSystem = " << dt.count() << "\n";
        //
        //--- Solve the linear system
        //
        //
        //--- Compute error norms
        //
        //
        //--- Output the approximation to GMSH
        //
        IMSI::OutputToGMSH("approximate.msh", myMesh, nullptr, myMesh.NumberVertices());
        //
    }
    Kokkos::finalize();
    return 0;
}
