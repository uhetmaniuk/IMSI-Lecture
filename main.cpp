//
// Created by Ulrich Hetmaniuk on 1/11/25.
//

#include <iostream>
#include <vector>

#include "src/Mesh.h"
#include "src/MeshUtils.h"

#include "Kokkos_Core.hpp"

int main(int argc, char *argv[]) {
    Kokkos::initialize(argc, argv);
    {

        IMSI::DomainParams dParams;
        dParams.numElePerDir[0] = 5;
        dParams.numElePerDir[1] = 4;

        auto myMesh = IMSI::GenerateMesh(dParams);
        IMSI::GetMatrixSparsity(myMesh);
    }
    Kokkos::finalize();
    return 0;
}