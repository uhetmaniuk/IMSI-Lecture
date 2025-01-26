#pragma once

#include <iostream>
#include <vector>

namespace IMSI {

    enum class RuleType : int {
        Gauss = 0, GaussLobatto
    };

    void getQuadrature
            (RuleType type,
             int sdim,
             int order,
             int &length,
             std::vector<double> &weights,
             std::vector<double> &xi,
             std::vector<double> &eta,
             std::vector<double> &zeta
            );

}