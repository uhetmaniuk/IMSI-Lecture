#pragma once

#include <functional>
#include <optional>
#include <tuple>

namespace IMSI {

    class Mesh;

    void MapDegreesOfFreedom(const std::vector<int>& bdyNodes,
                             std::vector<int>& globalToFree,
                             std::vector<int>& freeToGlobal);

//    std::tuple<double, double, double, double> GetErrorNorms( const Mesh& grid,
// double *u,
//                                                              const std::optional< std::function<double(double, double, double)> >& solution);

    void OutputToGMSH
            (
                    const char* fileName,
                    const Mesh& grid,
                    double *p,
                    int numDofs
            );

}
