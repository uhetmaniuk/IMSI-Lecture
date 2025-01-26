#pragma once

#include <array>

namespace IMSI {

  //! A class for two-dimensional Q2 finite element.
  
  /// The class defines the shape functions for the biquadratic isoparametric element.
  ///
  /// The nodes are numbered as follows
  ///
  ///    4 -7- 3
  ///
  ///    8  9  6
  ///
  ///    1 -5- 2
  ///
  /// We recommend a Gaussian quadrature with, at least, 2 points per direction.
  ///
  class fe2DQ2 {
    
  public:

      static const int sdim = 2;
      static const int numNode = 9;

      /// Constructor
    fe2DQ2() = default;

      /// \brief Returns the values of shape functions and their gradients at the quadrature point
      /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
      /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1]
      /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1] (unused)
      template<typename Scalar = double>
      static auto GetValuesGradients(Scalar xi, Scalar eta, [[maybe_unused]] Scalar zeta)  {
          std::array<Scalar, numNode * (sdim + 1)> values_gradients;
          int count = 0;
          //--------
          Scalar half(Scalar(0.5));
          Scalar one(Scalar(1.0));
          Scalar two(Scalar(2.0));
          //--------
          Scalar leftXi = half * (xi - one)*xi;
          Scalar middleXi = one - xi*xi;
          Scalar rightXi = half * (xi + one)*xi;
          //--------
          Scalar bottomEta = half * (eta - one)*eta;
          Scalar middleEta = one - eta*eta;
          Scalar topEta = half * (eta + one)*eta;
          //---------------------
          values_gradients[count++] = leftXi   * bottomEta;
          values_gradients[count++] = rightXi  * bottomEta;
          values_gradients[count++] = rightXi  * topEta;
          values_gradients[count++] = leftXi   * topEta;
          values_gradients[count++] = middleXi * bottomEta;
          values_gradients[count++] = rightXi  * middleEta;
          values_gradients[count++] = middleXi * topEta;
          values_gradients[count++] = leftXi   * middleEta;
          values_gradients[count++] = middleXi * middleEta;
          //
          values_gradients[count++] = (xi - half) * bottomEta;
          values_gradients[count++] = (xi + half) * bottomEta;
          values_gradients[count++] = (xi + half) * topEta;
          values_gradients[count++] = (xi - half) * topEta;
          values_gradients[count++] = - two * xi   * bottomEta;
          values_gradients[count++] = (xi + half) * middleEta;
          values_gradients[count++] = - two * xi   * topEta;
          values_gradients[count++] = (xi - half) * middleEta;
          values_gradients[count++] = - two * xi   * middleEta;
          //
          values_gradients[count++] = leftXi   * (eta - half);
          values_gradients[count++] = rightXi  * (eta - half);
          values_gradients[count++] = rightXi  * (eta + half);
          values_gradients[count++] = leftXi   * (eta + half);
          values_gradients[count++] = middleXi * (eta - half);
          values_gradients[count++] = rightXi  * (-two*eta);
          values_gradients[count++] = middleXi * (eta + half);
          values_gradients[count++] = leftXi   * (-two*eta);
          values_gradients[count++] = middleXi * (-two*eta);
          //
          return values_gradients;
      }

  };
  
}
