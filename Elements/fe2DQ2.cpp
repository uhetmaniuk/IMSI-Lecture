#include "Elements/fe2DQ2.h"

#include <cassert>
#include <cmath>
#include <iostream>


using namespace std;
using namespace FECode;


fe2DQ2::fe2DQ2
(
 const std::vector<double> &nodeCoord,
 const std::vector<double> &lQ,
 const std::vector<double> &wQ
)
: 
d_JxW(std::vector<double>(wQ.size())),
d_numGauss(wQ.size()),
d_Phi(std::vector<double>(9*wQ.size())),
d_PhiGrad(std::vector<double>(2*9*wQ.size())),
d_PointInt(std::vector<double>(2*wQ.size())),
d_x(std::vector<double>(9)),
d_y(std::vector<double>(9)),
d_lQuad(lQ),
d_wQuad(wQ),
d_hasShapeValue(false),
d_hasPointInt(false)
{
  
  size_t nPoint = 9;
  for (size_t j = 0; j < nPoint; ++j) 
  {
    d_x[j] = nodeCoord[2*j]; 
    d_y[j] = nodeCoord[2*j+1];
  }
  
  assert(d_numGauss >= 0);
  assert(d_lQuad.size() == 2 * d_wQuad.size());

}


void fe2DQ2::makeShapeValue() 
{
  
  size_t nPoint = 9;
  
  double xi, eta;
  double J[4], Jac = 0.0;
  double dqxi[9], dqeta[9];
  
  for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  {
    xi = d_lQuad[2*ipp];
    eta = d_lQuad[2*ipp + 1];
    //--------
    double leftXi = 0.5*(xi - 1.0)*xi;
    double middleXi = 1.0 - xi*xi;
    double rightXi = 0.5*(xi + 1.0)*xi;
    //--------
    double bottomEta = 0.5*(eta - 1.0)*eta;
    double middleEta = 1.0 - eta*eta;
    double topEta = 0.5*(eta + 1.0)*eta;
    //---------------------
    size_t pos = ipp;
    d_Phi[pos] = leftXi   * bottomEta; pos += d_numGauss;
    d_Phi[pos] = rightXi  * bottomEta; pos += d_numGauss;
    d_Phi[pos] = rightXi  * topEta;    pos += d_numGauss;
    d_Phi[pos] = leftXi   * topEta;    pos += d_numGauss;
    d_Phi[pos] = middleXi * bottomEta; pos += d_numGauss;
    d_Phi[pos] = rightXi  * middleEta; pos += d_numGauss;
    d_Phi[pos] = middleXi * topEta;    pos += d_numGauss;
    d_Phi[pos] = leftXi   * middleEta; pos += d_numGauss;
    d_Phi[pos] = middleXi * middleEta;
    //---------------------
    dqxi[0] = (xi - 0.5) * bottomEta;
    dqxi[1] = (xi + 0.5) * bottomEta;
    dqxi[2] = (xi + 0.5) * topEta;
    dqxi[3] = (xi - 0.5) * topEta;
    dqxi[4] = - 2.0*xi   * bottomEta;
    dqxi[5] = (xi + 0.5) * middleEta;
    dqxi[6] = - 2.0*xi   * topEta;
    dqxi[7] = (xi - 0.5) * middleEta;
    dqxi[8] = - 2.0*xi   * middleEta;
    //---------------------
    dqeta[0] = leftXi   * (eta - 0.5);
    dqeta[1] = rightXi  * (eta - 0.5);
    dqeta[2] = rightXi  * (eta + 0.5);
    dqeta[3] = leftXi   * (eta + 0.5);
    dqeta[4] = middleXi * (eta - 0.5);
    dqeta[5] = rightXi  * (-2.0*eta);
    dqeta[6] = middleXi * (eta + 0.5);
    dqeta[7] = leftXi   * (-2.0*eta);
    dqeta[8] = middleXi * (-2.0*eta);
    //---------------------
    J[0] = 0.0; J[1] = 0.0; J[2] = 0.0; J[3] = 0.0;
    for (size_t j = 0; j < nPoint; ++j) 
    {
      J[0] += d_x[j]*dqxi[j];
      J[1] += d_x[j]*dqeta[j];
      J[2] += d_y[j]*dqxi[j];
      J[3] += d_y[j]*dqeta[j];
    }
    //---------------------
    inverseJacobian(2, &(J[0]), Jac);
    if (Jac <= 0.0) 
    {
      if (Jac == 0.0)
        cerr << "\n fe2DQ2 >> Zero Jacobian determinant !!!\n\n";
      else
        cerr << "\n fe2DQ2 >> Negative Jacobian determinant !!!\n\n";
      bool hasPositiveJacobian = false;
      assert(hasPositiveJacobian == true);
    }
    //---------------------
    d_JxW[ipp] = d_wQuad[ipp] * Jac;
    //---------------------
    for (size_t j = 0; j < nPoint; ++j) {
      pos = 2*d_numGauss*j + 2*ipp;
      d_PhiGrad[pos]   = J[0]*dqxi[j] + J[2]*dqeta[j];
      d_PhiGrad[pos+1] = J[1]*dqxi[j] + J[3]*dqeta[j];
    }
    //---------------------
    
  } // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  
  d_hasShapeValue = true;
  
}


void fe2DQ2::makePointInt
(
) 
{
  
  if (d_hasShapeValue == false)
    makeShapeValue();
  
  size_t iN, iQ;
  for (iQ = 0; iQ < d_numGauss; ++iQ) 
  {
    double tmpX = 0.0, tmpY = 0.0;
    for (iN = 0; iN < 9; ++iN) {
      double shape = d_Phi[iN*d_numGauss + iQ];
      tmpX += d_x[iN]*shape;
      tmpY += d_y[iN]*shape;
    }
    unsigned int pos = 2*iQ;
    d_PointInt[pos] = tmpX;
    d_PointInt[pos+1] = tmpY;
  }
  
  d_hasPointInt = true;
  
}


void fe2DQ2::GetShapeGrad
(
 const size_t iP
 , const size_t iQ
 , std::valarray<double> &sg
)
{
  
  if (d_hasShapeValue == false)
    makeShapeValue();
  
  size_t pos = 2 * d_numGauss * iP + 2 * iQ;
  if (sg.size() != 2)
    sg.resize(2);
  sg[0] = d_PhiGrad[pos];
  sg[1] = d_PhiGrad[pos+1];
  
}


void fe2DQ2::GetShapeHessian
(
 const size_t iP
 , const size_t iQ
 , std::vector<double> &sh
)
{
  
  sh = std::vector<double>(9, 0.0);
  
  assert((iQ >= 0) && (iQ < d_numGauss));
  assert((iP >= 0) && (iP < 9));
  
  double xi = d_lQuad[2 * iQ];
  double eta = d_lQuad[2 * iQ + 1];
  
  double leftXi = 0.5*(xi - 1.0)*xi;
  double middleXi = 1.0 - xi*xi;
  double rightXi = 0.5*(xi + 1.0)*xi;
  //---------------------
  double bottomEta = 0.5*(eta - 1.0)*eta;
  double middleEta = 1.0 - eta*eta;
  double topEta = 0.5*(eta + 1.0)*eta;
  //---------------------
  double dqdxi[9];
  dqdxi[0] = (xi - 0.5) * bottomEta;
  dqdxi[1] = (xi + 0.5) * bottomEta;
  dqdxi[2] = (xi + 0.5) * topEta;
  dqdxi[3] = (xi - 0.5) * topEta;
  dqdxi[4] = - 2.0*xi   * bottomEta;
  dqdxi[5] = (xi + 0.5) * middleEta;
  dqdxi[6] = - 2.0*xi   * topEta;
  dqdxi[7] = (xi - 0.5) * middleEta;
  dqdxi[8] = - 2.0*xi   * middleEta;
  //---------------------
  double dqdeta[9];
  dqdeta[0] = leftXi   * (eta - 0.5);
  dqdeta[1] = rightXi  * (eta - 0.5);
  dqdeta[2] = rightXi  * (eta + 0.5);
  dqdeta[3] = leftXi   * (eta + 0.5);
  dqdeta[4] = middleXi * (eta - 0.5);
  dqdeta[5] = rightXi  * (-2.0*eta);
  dqdeta[6] = middleXi * (eta + 0.5);
  dqdeta[7] = leftXi   * (-2.0*eta);
  dqdeta[8] = middleXi * (-2.0*eta);
  //---------------------
  double d2qdxi2[9];
  d2qdxi2[0] = 1.0   * bottomEta;
  d2qdxi2[1] = 1.0   * bottomEta;
  d2qdxi2[2] = 1.0   * topEta;
  d2qdxi2[3] = 1.0   * topEta;
  d2qdxi2[4] = - 2.0 * bottomEta;
  d2qdxi2[5] = 1.0   * middleEta;
  d2qdxi2[6] = - 2.0 * topEta;
  d2qdxi2[7] = 1.0   * middleEta;
  d2qdxi2[8] = - 2.0 * middleEta;
  //---------------------
  double d2qdxideta[9];
  d2qdxideta[0] = (xi - 0.5) * (eta - 0.5);
  d2qdxideta[1] = (xi + 0.5) * (eta - 0.5);
  d2qdxideta[2] = (xi + 0.5) * (eta + 0.5);
  d2qdxideta[3] = (xi - 0.5) * (eta + 0.5);
  d2qdxideta[4] = (-2.0*xi)  * (eta - 0.5);
  d2qdxideta[5] = (xi + 0.5) * (-2.0*eta);
  d2qdxideta[6] = (-2.0*xi)  * (eta + 0.5);
  d2qdxideta[7] = (xi - 0.5) * (-2.0*eta);
  d2qdxideta[8] = (-2.0*xi)  * (-2.0*eta);
  //---------------------
  double d2qdeta2[9];
  d2qdeta2[0] = leftXi   * 1.0;
  d2qdeta2[1] = rightXi  * 1.0;
  d2qdeta2[2] = rightXi  * 1.0;
  d2qdeta2[3] = leftXi   * 1.0;
  d2qdeta2[4] = middleXi * 1.0;
  d2qdeta2[5] = rightXi  * (-2.0);
  d2qdeta2[6] = middleXi * 1.0;
  d2qdeta2[7] = leftXi   * (-2.0);
  d2qdeta2[8] = middleXi * (-2.0);
  //---------------------
  double dxdxi = 0.0, dydxi = 0.0, dxdeta = 0.0, dydeta = 0.0;
  double d2xdxi2 = 0.0, d2xdxideta = 0.0, d2xdeta2 = 0.0;
  double d2ydxi2 = 0.0, d2ydxideta = 0.0, d2ydeta2 = 0.0;
  for (unsigned int i = 0; i < 9; ++i) 
  {
    dxdxi  += d_x[i]*dqdxi[i];  dydxi  += d_y[i]*dqdxi[i];
    dxdeta += d_x[i]*dqdeta[i]; dydeta += d_y[i]*dqdeta[i];
    d2xdxi2 += d_x[i]*d2qdxi2[i]; d2ydxi2 += d_y[i]*d2qdxi2[i];
    d2xdxideta += d_x[i]*d2qdxideta[i]; d2ydxideta += d_y[i]*d2qdxideta[i];
    d2xdeta2 += d_x[i]*d2qdeta2[i]; d2ydeta2 += d_y[i]*d2qdeta2[i];
  }
  //---------------------
  double Jac = dxdxi * dydeta - dxdeta * dydxi;
  if (Jac <= 0.0) 
  {
    if (Jac == 0.0)
      cerr << "\n fe2DQ2::ShapeHessian >> Zero Jacobian determinant !!!\n\n";
    else
      cerr << "\n fe2DQ2::ShapeHessian >> Negative Jacobian determinant !!!\n\n";
    bool hasPositiveJacobian = false;
    assert(hasPositiveJacobian == true);
  }
  //---------------------
  Jac = 1.0 / Jac;
  double dqdx = (dydeta*dqdxi[iP] - dydxi*dqdeta[iP])*Jac;
  double dqdy = (-dxdeta*dqdxi[iP] + dxdxi*dqdeta[iP])*Jac;
  //---------------------
  d2qdxi2[iP] -= dqdx * d2xdxi2 + dqdy * d2ydxi2;
  d2qdxideta[iP] -= dqdx * d2xdxideta + dqdy * d2ydxideta;
  d2qdeta2[iP] -= dqdx * d2xdeta2 + dqdy * d2ydeta2;
  //---------------------
  double J11 = dxdxi*dxdxi, J12 = 2*dxdxi*dydxi, J13 = dydxi*dydxi;
  double J21 = dxdxi*dxdeta, J22 = dydxi*dxdeta + dxdxi*dydeta, J23 = dydxi*dydeta;
  double J31 = dxdeta*dxdeta, J32 = 2*dxdeta*dydeta, J33 = dydeta*dydeta;
  double iJ11 = J22*J33 - J32*J23;
  double iJ12 = J13*J32 - J33*J12;
  double iJ13 = J12*J23 - J22*J13;
  double iJ21 = J23*J31 - J33*J21;
  double iJ22 = J11*J33 - J31*J13;
  double iJ23 = J13*J21 - J23*J11;
  double iJ31 = J21*J32 - J31*J22;
  double iJ32 = J12*J31 - J32*J11;
  double iJ33 = J11*J22 - J21*J12;
  //---------------------
  Jac = J11*iJ11 + J12*iJ21 + J13*iJ31;
  //---------------------
  if (Jac == 0.0) 
  {
    cerr << "\n fe2DQ2::ShapeHessian >> Singular system !!!\n\n";
    bool isInvertibleForHessian = false;
    assert(isInvertibleForHessian == true);
  }
  //---------------------
  Jac = 1.0/Jac;
  sh[0] = Jac*(iJ11*d2qdxi2[iP] + iJ12*d2qdxideta[iP] + iJ13*d2qdeta2[iP]);
  sh[1] = Jac*(iJ21*d2qdxi2[iP] + iJ22*d2qdxideta[iP] + iJ23*d2qdeta2[iP]);
  sh[2] = sh[1];
  sh[3] = Jac*(iJ31*d2qdxi2[iP] + iJ32*d2qdxideta[iP] + iJ33*d2qdeta2[iP]);
  
}


double fe2DQ2::Volume
() 
{
  
  if (d_hasShapeValue == false)
    makeShapeValue();
  
  double vol = 0.0;
  for (int iG = 0; iG < d_numGauss; ++iG)
    vol += d_JxW[iG];
  
  return vol;
  
}


