#pragma once

#include <array>

namespace IMSI {
  
  //! A class for three-dimensional Q1 finite element.
  
  /// The class defines the shape functions for the trilinear isoparametric element.
  ///
  /// The nodes are numbered as follows
  ///
  ///    Front face
  ///
  ///    4 --- 3 
  ///
  ///    1 --- 2
  ///
  ///    Back face
  ///
  ///    8 --- 7 
  ///
  ///    5 --- 6
  ///
  /// We recommend a Gaussian quadrature with, at least, 2 points per direction.
  ///
  class fe3DQ1
  {
    
  public:

      static const int sdim = 3;
      static const int numNode = 8;

      /// Constructor
    fe3DQ1() = default;

      /// \brief Returns the values of shape functions and their gradients at the quadrature point
      /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
      /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1]
      /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1]
      template<typename Scalar>
      static auto GetValuesGradients(Scalar xi, Scalar eta, Scalar zeta)  {
          std::array<Scalar, numNode * (sdim + 1)> values_gradients;
          int count = 0;
          auto const one = Scalar(1.0);
          auto const oneOver8 = one / Scalar(8.0);
          //
          values_gradients[count++] = (one - xi) * (one - eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one - eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one + eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one + eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one - eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one - eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one + eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one + eta) * (one + zeta) * oneOver8;
          //
          values_gradients[count++] = - (one - eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one - eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one + eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = - (one + eta) * (one - zeta) * oneOver8;
          values_gradients[count++] = - (one - eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one - eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one + eta) * (one + zeta) * oneOver8;
          values_gradients[count++] = - (one + eta) * (one + zeta) * oneOver8;
          //
          values_gradients[count++] = - (one - xi) * (one - zeta) * oneOver8;
          values_gradients[count++] = - (one + xi) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one - zeta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one - zeta) * oneOver8;
          values_gradients[count++] = - (one - xi) * (one + zeta) * oneOver8;
          values_gradients[count++] = - (one + xi) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one + zeta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one + zeta) * oneOver8;
          //
          values_gradients[count++] = - (one - xi) * (one - eta) * oneOver8;
          values_gradients[count++] = - (one + xi) * (one - eta) * oneOver8;
          values_gradients[count++] = - (one + xi) * (one + eta) * oneOver8;
          values_gradients[count++] = - (one - xi) * (one + eta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one - eta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one - eta) * oneOver8;
          values_gradients[count++] = (one + xi) * (one + eta) * oneOver8;
          values_gradients[count++] = (one - xi) * (one + eta) * oneOver8;
          //
          return values_gradients;
      }

  };
  
}
