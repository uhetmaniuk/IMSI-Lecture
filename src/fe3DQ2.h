#pragma once

#include <array>

namespace IMSI {
  
  //! A class for three-dimensional Q2 finite element.
  
  /// The class defines the shape functions for the triquadratic isoparametric element.
  ///
  /// The node numbering is descriped on figure 3.7.6 in
  /// T. Hughes, "The Finite element method", Dover, 2000.
  ///
  /// We recommend a Gaussian quadrature with, at least, 3 points per direction.
  ///
  class fe3DQ2 {
    
  public:

      static const int sdim = 3;
      static const int numNode = 27;

      /// Constructor
      fe3DQ2() = default;

      /// \brief Returns the values of shape functions and their gradients at the quadrature point
      /// \param[in] xi 1st coordinate of quadrature point in [-1, 1]
      /// \param[in] eta 2nd coordinate of quadrature point in [-1, 1]
      /// \param[in] zeta 3rd coordinate of quadrature point in [-1, 1]
      template<typename Scalar>
      static auto GetValuesGradients(Scalar xi, Scalar eta, Scalar zeta) {
          std::array<Scalar, numNode * (sdim + 1)> values_gradients;
          int count = 0;
          auto const one = Scalar(1);
          auto const two = Scalar(2);
          auto const oneHalf = one / two;
          //---
          auto const leftXi = oneHalf*(xi - one)*xi;
          auto const middleXi = one - xi*xi;
          auto const rightXi = oneHalf*(xi + one)*xi;
          //---
          auto const bottomEta = oneHalf*(eta - one)*eta;
          auto const middleEta = one - eta*eta;
          auto const topEta = oneHalf*(eta + one)*eta;
          //---
          auto const downZeta = oneHalf*(zeta - one)*zeta;
          auto const middleZeta = one - zeta*zeta;
          auto const upZeta = oneHalf*(zeta + one)*zeta;
          //---------------------
          values_gradients[count++] = leftXi * bottomEta * downZeta;
          values_gradients[count++] = rightXi * bottomEta * downZeta;
          values_gradients[count++] = rightXi * topEta * downZeta;
          values_gradients[count++] = leftXi * topEta * downZeta;
          //
          values_gradients[count++] = leftXi * bottomEta * upZeta;
          values_gradients[count++] = rightXi * bottomEta * upZeta;
          values_gradients[count++] = rightXi * topEta * upZeta;
          values_gradients[count++] = leftXi * topEta * upZeta;
          //
          values_gradients[count++] = middleXi * bottomEta * downZeta;
          values_gradients[count++] = rightXi * middleEta * downZeta;
          values_gradients[count++] = middleXi * topEta * downZeta;
          values_gradients[count++] = leftXi * middleEta * downZeta;
          //
          values_gradients[count++] = middleXi * bottomEta * upZeta;
          values_gradients[count++] = rightXi * middleEta * upZeta;
          values_gradients[count++] = middleXi * topEta * upZeta;
          values_gradients[count++] = leftXi * middleEta * upZeta;
          //
          values_gradients[count++] = leftXi * bottomEta * middleZeta;
          values_gradients[count++] = rightXi * bottomEta * middleZeta;
          values_gradients[count++] = rightXi * topEta * middleZeta;
          values_gradients[count++] = leftXi * topEta * middleZeta;
          //
          values_gradients[count++] = middleXi * middleEta * downZeta;
          values_gradients[count++] = middleXi * middleEta * upZeta;
          //
          values_gradients[count++] = middleXi * bottomEta * middleZeta;
          values_gradients[count++] = middleXi * topEta * middleZeta;
          //
          values_gradients[count++] = leftXi * middleEta * middleZeta;
          values_gradients[count++] = rightXi * middleEta * middleZeta;
          //
          values_gradients[count++] = middleXi * middleEta * middleZeta;
          //
          //--- Derivative w.r.t. xi
          //
          values_gradients[count++] = (xi-oneHalf) * bottomEta * downZeta;
          values_gradients[count++] = (xi+oneHalf) * bottomEta * downZeta;
          values_gradients[count++] = (xi+oneHalf) * topEta * downZeta;
          values_gradients[count++] = (xi-oneHalf) * topEta * downZeta;
          //
          values_gradients[count++] = (xi-oneHalf) * bottomEta * upZeta;
          values_gradients[count++] = (xi+oneHalf) * bottomEta * upZeta;
          values_gradients[count++] = (xi+oneHalf) * topEta * upZeta;
          values_gradients[count++] = (xi-oneHalf) * topEta * upZeta;
          //
          values_gradients[count++] = (-two*xi) * bottomEta * downZeta;
          values_gradients[count++] = (xi+oneHalf) * middleEta * downZeta;
          values_gradients[count++] = (-two*xi) * topEta * downZeta;
          values_gradients[count++] = (xi-oneHalf) * middleEta * downZeta;
          //
          values_gradients[count++] = (-two*xi) * bottomEta * upZeta;
          values_gradients[count++] = (xi+oneHalf) * middleEta * upZeta;
          values_gradients[count++] = (-two*xi) * topEta * upZeta;
          values_gradients[count++] = (xi-oneHalf) * middleEta * upZeta;
          //
          values_gradients[count++] = (xi-oneHalf) * bottomEta * middleZeta;
          values_gradients[count++] = (xi+oneHalf) * bottomEta * middleZeta;
          values_gradients[count++] = (xi+oneHalf) * topEta * middleZeta;
          values_gradients[count++] = (xi-oneHalf) * topEta * middleZeta;
          //
          values_gradients[count++] = (-two*xi) * middleEta * downZeta;
          values_gradients[count++] = (-two*xi) * middleEta * upZeta;
          //
          values_gradients[count++] = (-two*xi) * bottomEta * middleZeta;
          values_gradients[count++] = (-two*xi) * topEta * middleZeta;
          //
          values_gradients[count++] = (xi-oneHalf) * middleEta * middleZeta;
          values_gradients[count++] = (xi+oneHalf) * middleEta * middleZeta;
          //
          values_gradients[count++] = (-two*xi) * middleEta * middleZeta;
          //
          //--- Derivative w.r.t. eta
          //
          values_gradients[count++] = leftXi * (eta-oneHalf) * downZeta;
          values_gradients[count++] = rightXi * (eta-oneHalf) * downZeta;
          values_gradients[count++] = rightXi * (eta+oneHalf) * downZeta;
          values_gradients[count++] = leftXi * (eta+oneHalf) * downZeta;
          //
          values_gradients[count++] = leftXi * (eta-oneHalf) * upZeta;
          values_gradients[count++] = rightXi * (eta-oneHalf) * upZeta;
          values_gradients[count++] = rightXi * (eta+oneHalf) * upZeta;
          values_gradients[count++] = leftXi * (eta+oneHalf) * upZeta;
          //
          values_gradients[count++] = middleXi * (eta-oneHalf) * downZeta;
          values_gradients[count++] = rightXi * (-two*eta) * downZeta;
          values_gradients[count++] = middleXi * (eta+oneHalf) * downZeta;
          values_gradients[count++] = leftXi * (-two*eta) * downZeta;
          //
          values_gradients[count++] = middleXi * (eta-oneHalf) * upZeta;
          values_gradients[count++] = rightXi * (-two*eta) * upZeta;
          values_gradients[count++] = middleXi * (eta+oneHalf) * upZeta;
          values_gradients[count++] = leftXi * (-two*eta) * upZeta;
          //
          values_gradients[count++] = leftXi * (eta-oneHalf) * middleZeta;
          values_gradients[count++] = rightXi * (eta-oneHalf) * middleZeta;
          values_gradients[count++] = rightXi * (eta+oneHalf) * middleZeta;
          values_gradients[count++] = leftXi * (eta+oneHalf) * middleZeta;
          //
          values_gradients[count++] = middleXi * (-two*eta) * downZeta;
          values_gradients[count++] = middleXi * (-two*eta) * upZeta;
          //
          values_gradients[count++] = middleXi * (eta-oneHalf) * middleZeta;
          values_gradients[count++] = middleXi * (eta+oneHalf) * middleZeta;
          //
          values_gradients[count++] = leftXi * (-two*eta) * middleZeta;
          values_gradients[count++] = rightXi * (-two*eta) * middleZeta;
          //
          values_gradients[count++] = middleXi * (-two*eta) * middleZeta;
          //
          //--- Derivative w.r.t. zeta
          //
          values_gradients[count++] = leftXi * bottomEta * (zeta-oneHalf);
          values_gradients[count++] = rightXi * bottomEta * (zeta-oneHalf);
          values_gradients[count++] = rightXi * topEta * (zeta-oneHalf);
          values_gradients[count++] = leftXi * topEta * (zeta-oneHalf);
          //
          values_gradients[count++] = leftXi * bottomEta * (zeta+oneHalf);
          values_gradients[count++] = rightXi * bottomEta * (zeta+oneHalf);
          values_gradients[count++] = rightXi * topEta * (zeta+oneHalf);
          values_gradients[count++] = leftXi * topEta * (zeta+oneHalf);
          //
          values_gradients[count++] = middleXi * bottomEta * (zeta-oneHalf);
          values_gradients[count++] = rightXi * middleEta * (zeta-oneHalf);
          values_gradients[count++] = middleXi * topEta * (zeta-oneHalf);
          values_gradients[count++] = leftXi * middleEta * (zeta-oneHalf);
          //
          values_gradients[count++] = middleXi * bottomEta * (zeta+oneHalf);
          values_gradients[count++] = rightXi * middleEta * (zeta+oneHalf);
          values_gradients[count++] = middleXi * topEta * (zeta+oneHalf);
          values_gradients[count++] = leftXi * middleEta * (zeta+oneHalf);
          //
          values_gradients[count++] = leftXi * bottomEta * (-two*zeta);
          values_gradients[count++] = rightXi * bottomEta * (-two*zeta);
          values_gradients[count++] = rightXi * topEta * (-two*zeta);
          values_gradients[count++] = leftXi * topEta * (-two*zeta);
          //
          values_gradients[count++] = middleXi * middleEta * (zeta-oneHalf);
          values_gradients[count++] = middleXi * middleEta * (zeta+oneHalf);
          //
          values_gradients[count++] = middleXi * bottomEta * (-two*zeta);
          values_gradients[count++] = middleXi * topEta * (-two*zeta);
          //
          values_gradients[count++] = leftXi * middleEta * (-two*zeta);
          values_gradients[count++] = rightXi * middleEta * (-two*zeta);
          //
          values_gradients[count++] = middleXi * middleEta * (-two*zeta);
          //
          return values_gradients;
      }

      };
  
}
