#pragma once

#include <Kokkos_Core.hpp>

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
  KOKKOS_INLINE_FUNCTION
  fe1DQ1() = default;

  /// \brief Returns the values of shape functions and their gradients at the quadrature point
  /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
  /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1] (unused)
  /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  static auto
  GetValuesGradients(
      Scalar                  xi,
      [[maybe_unused]] Scalar eta,
      [[maybe_unused]] Scalar zeta,
      Scalar*                 values_gradients)
  {
    values_gradients[0] = (1. - xi) * 0.5;
    values_gradients[1] = (1. + xi) * 0.5;
    //
    values_gradients[2] = -0.5;
    values_gradients[3] = 0.5;
  }
};

}  // namespace IMSI
