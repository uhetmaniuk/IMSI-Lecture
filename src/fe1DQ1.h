#pragma once

#include <array>

namespace IMSI {

//! A class for one-dimensional P1 finite element.

/// The class defines the shape functions for the linear isoparametric element.
///
/// The nodes are numbered as follows
///
///    1 -- 2
///
/// We recommend a Gaussian quadrature with, at least, 2 points per direction.
///
class fe1DQ1
{
 public:
  static const int sdim    = 1;
  static const int numNode = 2;

  /// Constructor
  fe1DQ1() = default;

  /// \brief Returns the values of shape functions and their gradients at the quadrature point
  /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
  /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1] (unused)
  /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
  static auto
  GetValuesGradients(double xi, [[maybe_unused]] double eta, [[maybe_unused]] double zeta)
  {
    std::array<double, numNode*(sdim + 1)> values_gradients;
    int                                    count = 0;
    values_gradients[count++]                    = (1. - xi) * 0.5;
    values_gradients[count++]                    = (1. + xi) * 0.5;
    //
    values_gradients[count++] = -0.5;
    values_gradients[count++] = 0.5;
    //
    return values_gradients;
  }
};

}  // namespace IMSI
