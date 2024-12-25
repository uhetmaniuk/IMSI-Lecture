#include "Elements/fe3DQ2.h"

#include <cassert>
#include <cmath>
#include <iostream>


using namespace std;
using namespace FECode;


fe3DQ2::fe3DQ2
(
 const std::vector<double> &nodeCoord,
 const std::vector<double> &lQ,
 const std::vector<double> &wQ
)
: 
d_JxW(vector<double>(wQ.size())),
d_numGauss(wQ.size()),
d_Phi(std::vector<double>(27*wQ.size())),
d_PhiGrad(std::vector<double>(3*27*wQ.size())),
d_PointInt(std::vector<double>(3*wQ.size())),
d_lQuad(lQ),
d_wQuad(wQ)
{

  size_t nPoint = 27;
  double x[27], y[27], z[27];
  for (unsigned int j = 0; j < nPoint; ++j) {
    x[j] = nodeCoord[3*j];
    y[j] = nodeCoord[3*j+1];
    z[j] = nodeCoord[3*j+2];
  }

  double xi, eta, zeta;
  double J[9], Jac = 0.0;
  double dqxi[27], dqeta[27], dqzeta[27];

  for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  {
    xi = d_lQuad[3*ipp];
    eta = d_lQuad[3*ipp + 1];
    zeta = d_lQuad[3*ipp + 2];
    //---
    double leftXi = 0.5*(xi - 1.0)*xi;
    double middleXi = 1.0 - xi*xi;
    double rightXi = 0.5*(xi + 1.0)*xi;
    //---
      double bottomEta = 0.5*(eta - 1.0)*eta;
      double middleEta = 1.0 - eta*eta;
      double topEta = 0.5*(eta + 1.0)*eta;
      //---
        double downZeta = 0.5*(zeta - 1.0)*zeta;
        double middleZeta = 1.0 - zeta*zeta;
        double upZeta = 0.5*(zeta + 1.0)*zeta;
        //---------------------
        size_t pos = ipp;
        d_Phi[pos] = leftXi * bottomEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * bottomEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * topEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = leftXi * topEta * downZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = leftXi * bottomEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * bottomEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * topEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = leftXi * topEta * upZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = middleXi * bottomEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * middleEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = middleXi * topEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = leftXi * middleEta * downZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = middleXi * bottomEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * middleEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = middleXi * topEta * upZeta;
        pos += d_numGauss; d_Phi[pos] = leftXi * middleEta * upZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = leftXi * bottomEta * middleZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * bottomEta * middleZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * topEta * middleZeta;
        pos += d_numGauss; d_Phi[pos] = leftXi * topEta * middleZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = middleXi * middleEta * downZeta;
        pos += d_numGauss; d_Phi[pos] = middleXi * middleEta * upZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = middleXi * bottomEta * middleZeta;
        pos += d_numGauss; d_Phi[pos] = middleXi * topEta * middleZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = leftXi * middleEta * middleZeta;
        pos += d_numGauss; d_Phi[pos] = rightXi * middleEta * middleZeta;
        //---
        pos += d_numGauss; d_Phi[pos] = middleXi * middleEta * middleZeta;
        //---------------------
        dqxi[ 0] = (xi-0.5) * bottomEta * downZeta;
        dqxi[ 1] = (xi+0.5) * bottomEta * downZeta;
        dqxi[ 2] = (xi+0.5) * topEta * downZeta;
        dqxi[ 3] = (xi-0.5) * topEta * downZeta;
        //---
        dqxi[ 4] = (xi-0.5) * bottomEta * upZeta;
        dqxi[ 5] = (xi+0.5) * bottomEta * upZeta;
        dqxi[ 6] = (xi+0.5) * topEta * upZeta;
        dqxi[ 7] = (xi-0.5) * topEta * upZeta;
        //---
        dqxi[ 8] = (-2.0*xi) * bottomEta * downZeta;
        dqxi[ 9] = (xi+0.5) * middleEta * downZeta;
        dqxi[10] = (-2.0*xi) * topEta * downZeta;
        dqxi[11] = (xi-0.5) * middleEta * downZeta;
        //---
        dqxi[12] = (-2.0*xi) * bottomEta * upZeta;
        dqxi[13] = (xi+0.5) * middleEta * upZeta;
        dqxi[14] = (-2.0*xi) * topEta * upZeta;
        dqxi[15] = (xi-0.5) * middleEta * upZeta;
        //---
        dqxi[16] = (xi-0.5) * bottomEta * middleZeta;
        dqxi[17] = (xi+0.5) * bottomEta * middleZeta;
        dqxi[18] = (xi+0.5) * topEta * middleZeta;
        dqxi[19] = (xi-0.5) * topEta * middleZeta;
        //---
        dqxi[20] = (-2.0*xi) * middleEta * downZeta;
        dqxi[21] = (-2.0*xi) * middleEta * upZeta;
        //---
        dqxi[22] = (-2.0*xi) * bottomEta * middleZeta;
        dqxi[23] = (-2.0*xi) * topEta * middleZeta;
        //---
        dqxi[24] = (xi-0.5) * middleEta * middleZeta;
        dqxi[25] = (xi+0.5) * middleEta * middleZeta;
        //---
        dqxi[26] = (-2.0*xi) * middleEta * middleZeta;
        //---------------------
        dqeta[ 0] = leftXi * (eta-0.5) * downZeta;
        dqeta[ 1] = rightXi * (eta-0.5) * downZeta;
        dqeta[ 2] = rightXi * (eta+0.5) * downZeta;
        dqeta[ 3] = leftXi * (eta+0.5) * downZeta;
        //---
        dqeta[ 4] = leftXi * (eta-0.5) * upZeta;
        dqeta[ 5] = rightXi * (eta-0.5) * upZeta;
        dqeta[ 6] = rightXi * (eta+0.5) * upZeta;
        dqeta[ 7] = leftXi * (eta+0.5) * upZeta;
        //---
        dqeta[ 8] = middleXi * (eta-0.5) * downZeta;
        dqeta[ 9] = rightXi * (-2.0*eta) * downZeta;
        dqeta[10] = middleXi * (eta+0.5) * downZeta;
        dqeta[11] = leftXi * (-2.0*eta) * downZeta;
        //---
        dqeta[12] = middleXi * (eta-0.5) * upZeta;
        dqeta[13] = rightXi * (-2.0*eta) * upZeta;
        dqeta[14] = middleXi * (eta+0.5) * upZeta;
        dqeta[15] = leftXi * (-2.0*eta) * upZeta;
        //---
        dqeta[16] = leftXi * (eta-0.5) * middleZeta;
        dqeta[17] = rightXi * (eta-0.5) * middleZeta;
        dqeta[18] = rightXi * (eta+0.5) * middleZeta;
        dqeta[19] = leftXi * (eta+0.5) * middleZeta;
        //---
        dqeta[20] = middleXi * (-2.0*eta) * downZeta;
        dqeta[21] = middleXi * (-2.0*eta) * upZeta;
        //---
        dqeta[22] = middleXi * (eta-0.5) * middleZeta;
        dqeta[23] = middleXi * (eta+0.5) * middleZeta;
        //---
        dqeta[24] = leftXi * (-2.0*eta) * middleZeta;
        dqeta[25] = rightXi * (-2.0*eta) * middleZeta;
        //---
        dqeta[26] = middleXi * (-2.0*eta) * middleZeta;
        //---------------------
        dqzeta[ 0] = leftXi * bottomEta * (zeta-0.5);
        dqzeta[ 1] = rightXi * bottomEta * (zeta-0.5);
        dqzeta[ 2] = rightXi * topEta * (zeta-0.5);
        dqzeta[ 3] = leftXi * topEta * (zeta-0.5);
        //---
        dqzeta[ 4] = leftXi * bottomEta * (zeta+0.5);
        dqzeta[ 5] = rightXi * bottomEta * (zeta+0.5);
        dqzeta[ 6] = rightXi * topEta * (zeta+0.5);
        dqzeta[ 7] = leftXi * topEta * (zeta+0.5);
        //---
        dqzeta[ 8] = middleXi * bottomEta * (zeta-0.5);
        dqzeta[ 9] = rightXi * middleEta * (zeta-0.5);
        dqzeta[10] = middleXi * topEta * (zeta-0.5);
        dqzeta[11] = leftXi * middleEta * (zeta-0.5);
        //---
        dqzeta[12] = middleXi * bottomEta * (zeta+0.5);
        dqzeta[13] = rightXi * middleEta * (zeta+0.5);
        dqzeta[14] = middleXi * topEta * (zeta+0.5);
        dqzeta[15] = leftXi * middleEta * (zeta+0.5);
        //---
        dqzeta[16] = leftXi * bottomEta * (-2.0*zeta);
        dqzeta[17] = rightXi * bottomEta * (-2.0*zeta);
        dqzeta[18] = rightXi * topEta * (-2.0*zeta);
        dqzeta[19] = leftXi * topEta * (-2.0*zeta);
        //---
        dqzeta[20] = middleXi * middleEta * (zeta-0.5);
        dqzeta[21] = middleXi * middleEta * (zeta+0.5);
        //---
        dqzeta[22] = middleXi * bottomEta * (-2.0*zeta);
        dqzeta[23] = middleXi * topEta * (-2.0*zeta);
        //---
        dqzeta[24] = leftXi * middleEta * (-2.0*zeta);
        dqzeta[25] = rightXi * middleEta * (-2.0*zeta);
        //---
        dqzeta[26] = middleXi * middleEta * (-2.0*zeta);
        //---------------------
        J[0] = 0.0; J[3] = 0.0; J[6] = 0.0;
        J[1] = 0.0; J[4] = 0.0; J[7] = 0.0;
        J[2] = 0.0; J[5] = 0.0; J[8] = 0.0;
        for (size_t j = 0; j < nPoint; ++j){
          J[0] += x[j]*dqxi[j]; J[3] += y[j]*dqxi[j]; J[6] += z[j]*dqxi[j];
          J[1] += x[j]*dqeta[j]; J[4] += y[j]*dqeta[j]; J[7] += z[j]*dqeta[j];
          J[2] += x[j]*dqzeta[j]; J[5] += y[j]*dqzeta[j]; J[8] += z[j]*dqzeta[j];
        }
        //---------------------
        inverseJacobian(3, &(J[0]), Jac);
        if (Jac <= 0.0) {
          if (Jac == 0.0)
            cerr << "\n fe3DQ2 >> Zero Jacobian determinant !!!\n\n";
          else
            cerr << "\n fe3DQ2 >> Negative Jacobian determinant !!!\n\n";
          bool hasPositiveJacobian = false;
          assert(hasPositiveJacobian == true);
        }
        //---------------------
        d_JxW[ipp] = d_wQuad[ipp] * Jac;
        //---------------------
        for (size_t j = 0; j < nPoint; ++j) {
          pos = 3 * d_numGauss * j + 3 * ipp;
          d_PhiGrad[pos]   = J[0]*dqxi[j] + J[3]*dqeta[j] + J[6]*dqzeta[j];
          d_PhiGrad[pos+1] = J[1]*dqxi[j] + J[4]*dqeta[j] + J[7]*dqzeta[j];
          d_PhiGrad[pos+2] = J[2]*dqxi[j] + J[5]*dqeta[j] + J[8]*dqzeta[j];
        }
        //---------------------
        double tmpX = 0.0, tmpY = 0.0, tmpZ = 0.0;
        for (size_t j = 0; j < nPoint; ++j) {
          double shape = d_Phi[j*d_numGauss + ipp];
          tmpX += x[j]*shape;
          tmpY += y[j]*shape;
          tmpZ += z[j]*shape;
        }
        pos = 3 * ipp;
        d_PointInt[pos] = tmpX;
        d_PointInt[pos+1] = tmpY;
        d_PointInt[pos+2] = tmpZ;
        //---------------------
    
  } // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

}


void fe3DQ2::GetShapeGrad
(
  const size_t iP
  , const size_t iQ
  , std::valarray<double> &sg
)
{

  size_t pos = 3 * d_numGauss * iP + 3 * iQ;
  if (sg.size() != 3)
    sg.resize(3);
  sg[0] = d_PhiGrad[pos];
  sg[1] = d_PhiGrad[pos + 1];
  sg[2] = d_PhiGrad[pos + 2];

}


void fe3DQ2::GetShapeHessian
(
  const size_t iP
  , const size_t iQ
  , std::vector<double> &sh
)
{

  sh = std::vector<double>(27, 0.0);

}


double fe3DQ2::Volume() {

  double vol = 0.0;
  for (int iG = 0; iG < d_numGauss; ++iG)
    vol += d_JxW[iG];

  return vol;

}


