#pragma once

#include <array>

namespace IMSI {
  
  //! A class for two-dimensional Q1 finite element.
  
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
  class fe2DQ1 {
    
  public:

      static const int sdim = 2;
      static const int numNode = 4;

    /// Constructor
    fe2DQ1() = default;

    /// \brief Returns the values of shape functions and their gradients at the quadrature point
    /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
    /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1]
    /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
    template<typename Scalar>
	static auto GetValuesGradients(Scalar xi, Scalar eta, [[maybe_unused]] Scalar zeta)  {
        std::array<Scalar, numNode * (sdim + 1)> values_gradients;
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
        return values_gradients;
    }

  };
  
}
