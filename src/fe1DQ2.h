#pragma once

#include <array>

namespace IMSI {
  
  //! A class for one-dimensional Q2 finite element.
  
  /// The class defines the shape functions for the quadratic isoparametric element.
  ///
  /// The nodes are numbered as follows
  ///
  ///    1 -3- 2
  ///
  /// We recommend a Gaussian quadrature with, at least, 3 points per direction.
  ///
  class fe1DQ2 {

  public:

      static const int sdim = 1;
      static const int numNode = 3;

      /// Constructor
      fe1DQ2() = default;

      /// \brief Returns the values of shape functions and their gradients at the quadrature point
      /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
      /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1] (unused)
      /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
      static auto GetValuesGradients(double xi, double eta, [[maybe_unused]] double zeta)  {
          std::array<double, numNode * (sdim + 1)> values_gradients;
          int count = 0;
          values_gradients[count++] = 0.5*(xi - 1.0)*xi;
          values_gradients[count++] = 0.5*(xi + 1.0)*xi;
          values_gradients[count++] = 1.0 - xi*xi;
          //
          values_gradients[count++] = (xi - 0.5);
          values_gradients[count++] = (xi + 0.5);
          values_gradients[count++] = - 2.0 * xi;
          return values_gradients;
      }

  };
  
}
