#pragma once

#include <functional>
#include <optional>

namespace IMSI {

    const std::function<double(double, double, double)> constantOne = [](double x, double y, double) -> double {
        return double(1.0);
    };

    /// \brief Set of functions for 2D basic Poisson problem
    struct ParabolicPb {
        const std::function<double(double, double, double)> f= [](double x, double y, double) -> double {
            return double(2.0 * x * (1.0 - x) + 2.0 * y * (1.0 - y)); };
        const std::function<double(double, double, double)> ax = constantOne;
        const std::function<double(double, double, double)> ay = constantOne;
        const std::optional< std::function<double(double, double, double)> > u = [](double x, double y, double) -> double {
            return double(x * (1.0 - x) * y * (1.0 - y)); };
    };

}
