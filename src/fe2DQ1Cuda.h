#pragma once

#include <Kokkos_Core.hpp>

namespace IMSI {

  //! A class for two-dimensional Q1 finite element (CUDA-compatible)

  /// The class defines the shape functions for the bilinear isoparametric element.
  ///
  /// The nodes are numbered as follows
  ///
  ///    4 --- 3
  ///
  ///    1 --- 2
  ///
  /// We recommend a Gaussian quadrature with, at least, 2 points per direction.
  ///
  class fe2DQ1Cuda {

  public:

      static constexpr int sdim = 2;
      static constexpr int numNode = 4;

    /// Constructor
    KOKKOS_INLINE_FUNCTION
    fe2DQ1Cuda() = default;

    /// \brief Returns the values of shape functions and their gradients at the quadrature point
    /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
    /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1]
    /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
    /// \param[out] values_gradients Array of size numNode * (sdim + 1) to store results
    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    static void GetValuesGradients(Scalar xi, Scalar eta, [[maybe_unused]] Scalar zeta,
                                   Scalar* values_gradients)  {
        int count = 0;
        auto const one = Scalar(1.0);
        auto const oneOver4 = one / Scalar(4.0);
        //
        values_gradients[count++] = (one - xi) * (one - eta) * oneOver4;
        values_gradients[count++] = (one + xi) * (one - eta) * oneOver4;
        values_gradients[count++] = (one + xi) * (one + eta) * oneOver4;
        values_gradients[count++] = (one - xi) * (one + eta) * oneOver4;
        //
        values_gradients[count++] = - (one - eta) * oneOver4;
        values_gradients[count++] = (one - eta) * oneOver4;
        values_gradients[count++] = (one + eta) * oneOver4;
        values_gradients[count++] = - (one + eta) * oneOver4;
        //
        values_gradients[count++] = - (one - xi) * oneOver4;
        values_gradients[count++] = - (one + xi) * oneOver4;
        values_gradients[count++] = (one + xi) * oneOver4;
        values_gradients[count++] = (one - xi) * oneOver4;
    }

  };

}
