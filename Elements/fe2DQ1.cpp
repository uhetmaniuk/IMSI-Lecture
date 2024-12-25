#include "Elements/fe2DQ1.h"

#include <cassert>
#include <cmath>
#include <iostream>


using namespace std;
using namespace FECode;


fe2DQ1::fe2DQ1
(
 const std::vector<double> &nodeCoord,
 const std::vector<double> &lQ,
 const std::vector<double> &wQ
)
: 
d_JxW(vector<double>(wQ.size())),
d_numGauss(wQ.size()),
d_Phi(std::vector<double>(4*wQ.size())),
d_PhiGrad(std::vector<double>(2*4*wQ.size())),
d_PointInt(std::vector<double>(2*wQ.size())),
d_x(std::vector<double>(4)),
d_y(std::vector<double>(4)),
d_lQuad(lQ),
d_wQuad(wQ),
d_hasShapeValue(false),
d_hasPointInt(false)
{
  
  d_x[0] = nodeCoord[0]; d_y[0] = nodeCoord[1];
  d_x[1] = nodeCoord[2]; d_y[1] = nodeCoord[3];
  d_x[2] = nodeCoord[4]; d_y[2] = nodeCoord[5];
  d_x[3] = nodeCoord[6]; d_y[3] = nodeCoord[7];
  
}


void fe2DQ1::makeShapeValue
(
) 
{
  
  double xi, eta;
  double J[4], Jac = 0.0;
  
  double dqxi[4], dqeta[4];
  
  for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  {
    xi = d_lQuad[2*ipp];
    eta = d_lQuad[2*ipp + 1];
    //--------
    double halfOnePXi = 0.5*(1.0 + xi);
    double halfOneMXi = 1.0 - halfOnePXi;
    //--------
    double halfOnePEta = 0.5*(1.0 + eta);
    double halfOneMEta = 1.0 - halfOnePEta;
    //---------------------
    unsigned int pos = ipp;
    d_Phi[pos] = halfOneMXi * halfOneMEta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOneMEta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOnePEta; pos += d_numGauss;
    d_Phi[pos] = halfOneMXi * halfOnePEta;
    //---------------------
    dqxi[0] = -0.5*halfOneMEta;
    dqxi[1] = 0.5*halfOneMEta;
    dqxi[2] = 0.5*halfOnePEta;
    dqxi[3] = -0.5*halfOnePEta;
    //---------------------
    dqeta[0] = -0.5*halfOneMXi;
    dqeta[1] = -0.5*halfOnePXi;
    dqeta[2] = 0.5*halfOnePXi;
    dqeta[3] = 0.5*halfOneMXi;
    //---------------------
    // Define the matrix J by writing
    //
    // [1, 0; 0, 1] = [dxi/dx, deta/dx; dxi/dy, deta/dy]
    //                [dx/dxi, dy/dxi; dx/deta, dy/deta]
    //
    // and J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
    //
    // x = \sum_{A} x_{A} q_{A}(xi, eta)
    // y = \sum_{A} y_{A} q_{A}(xi, eta)
    //
    J[0] = d_x[0]*dqxi[0] + d_x[1]*dqxi[1] + d_x[2]*dqxi[2] + d_x[3]*dqxi[3];
    J[1] = d_x[0]*dqeta[0] + d_x[1]*dqeta[1] + d_x[2]*dqeta[2] + d_x[3]*dqeta[3];
    J[2] = d_y[0]*dqxi[0] + d_y[1]*dqxi[1] + d_y[2]*dqxi[2] + d_y[3]*dqxi[3];
    J[3] = d_y[0]*dqeta[0] + d_y[1]*dqeta[1] + d_y[2]*dqeta[2] + d_y[3]*dqeta[3];
    //---------------------
    inverseJacobian(2, &(J[0]), Jac);
    if (Jac <= 0.0) {
      if (Jac == 0.0)
        cerr << "\n fe2DQ1 >> Zero Jacobian determinant !!!\n\n";
      else
        cerr << "\n fe2DQ1 >> Negative Jacobian determinant !!!\n\n";
      bool hasPositiveJacobian = false;
      assert(hasPositiveJacobian == true);
    }
    //---------------------
    d_JxW[ipp] = d_wQuad[ipp] * Jac;
    //---------------------
    // Do the operations
    //
    // q(xi, eta) = phi(x(xi,eta), y(xi,eta))
    //
    // dq/dxi = dphi/dx dx/dxi + dphi/dy dy/dxi
    // dq/deta = dphi/dx dx/deta + dphi/dy dy/deta
    //
    // [dq/dxi; dq/deta] = J [dphi/dx; dphi/dy]
    // (where J is defined as [dx/dxi, dy/dxi; dx/deta, dy/deta])
    //
    // The variable J contains "J^{-1}" after 'inverseJacobian'
    //
    for (unsigned int j = 0; j < nnode; ++j) {
      pos = sdim*d_numGauss*j + sdim*ipp;
      d_PhiGrad[pos]   = J[0]*dqxi[j] + J[2]*dqeta[j];
      d_PhiGrad[pos+1] = J[1]*dqxi[j] + J[3]*dqeta[j];
    }
    //---------------------
    
  } // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  
  d_hasShapeValue = true;
  
}


void fe2DQ1::makePointInt
(
) 
{
  
  if (d_hasShapeValue == false)
    makeShapeValue();
  
  unsigned int iN, iQ;
  for (iQ = 0; iQ < d_numGauss; ++iQ) {
    double tmpX = 0.0, tmpY = 0.0;
    for (iN = 0; iN < 4; ++iN) {
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



void fe2DQ1::GetShapeGrad
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


void fe2DQ1::GetShapeHessian
(
 const size_t iP
 , const size_t iQ
 , std::vector<double> &sh
)
{
  
  sh = std::vector<double>(4, 0.0);
  
  if ((iQ < 0) || (iQ >= d_numGauss))
    return;
  
  if ((iP < 0) || (iP >= 4))
    return;
  
  double xi = d_lQuad[2 * iQ];
  double eta = d_lQuad[2 * iQ + 1];
  
  double halfOnePXi = 0.5*(1.0 + xi);
  double halfOneMXi = 1.0 - halfOnePXi;
  double halfOnePEta = 0.5*(1.0 + eta);
  double halfOneMEta = 1.0 - halfOnePEta;
  //---------------------
  double dqdxi[4];
  dqdxi[0] = (-0.5) * halfOneMEta;
  dqdxi[1] = ( 0.5) * halfOneMEta;
  dqdxi[2] = ( 0.5) * halfOnePEta;
  dqdxi[3] = (-0.5) * halfOnePEta;
  //---------------------
  double dqdeta[4];
  dqdeta[0] = halfOneMXi * (-0.5);
  dqdeta[1] = halfOnePXi * (-0.5);
  dqdeta[2] = halfOnePXi * (0.5);
  dqdeta[3] = halfOneMXi * (0.5);
  //---------------------
  double d2qdxi2[4];
  d2qdxi2[0] = 0.0;
  d2qdxi2[1] = 0.0;
  d2qdxi2[2] = 0.0;
  d2qdxi2[3] = 0.0;
  //---------------------
  double d2qdxideta[4];
  d2qdxideta[0] = 0.25;
  d2qdxideta[1] = -0.25;
  d2qdxideta[2] = 0.25;
  d2qdxideta[3] = -0.25;
  //---------------------
  double d2qdeta2[4];
  d2qdeta2[0] = 0.0;
  d2qdeta2[1] = 0.0;
  d2qdeta2[2] = 0.0;
  d2qdeta2[3] = 0.0;
  //---------------------
  double dxdxi = 0.0, dydxi = 0.0, dxdeta = 0.0, dydeta = 0.0;
  double d2xdxi2 = 0.0, d2xdxideta = 0.0, d2xdeta2 = 0.0;
  double d2ydxi2 = 0.0, d2ydxideta = 0.0, d2ydeta2 = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    dxdxi  += d_x[i]*dqdxi[i];  dydxi  += d_y[i]*dqdxi[i];
    dxdeta += d_x[i]*dqdeta[i]; dydeta += d_y[i]*dqdeta[i];
    d2xdxi2 += d_x[i]*d2qdxi2[i]; d2ydxi2 += d_y[i]*d2qdxi2[i];
    d2xdxideta += d_x[i]*d2qdxideta[i]; d2ydxideta += d_y[i]*d2qdxideta[i];
    d2xdeta2 += d_x[i]*d2qdeta2[i]; d2ydeta2 += d_y[i]*d2qdeta2[i];
  }
  //---------------------
  double Jac = dxdxi * dydeta - dxdeta * dydxi;
  if (Jac <= 0.0) {
    if (Jac == 0.0)
      cerr << "\n fe2DQ1::ShapeHessian >> Zero Jacobian determinant !!!\n\n";
    else
      cerr << "\n fe2DQ1::ShapeHessian >> Negative Jacobian determinant !!!\n\n";
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
  if (Jac == 0.0) {
    cerr << "\n fe2DQ1::ShapeHessian >> Singular system !!!\n\n";
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


double fe2DQ1::Volume
(
) 
{
  
  if (d_hasShapeValue == false)
    makeShapeValue();
  
  double vol = 0.0;
  for (int iG = 0; iG < d_numGauss; ++iG)
    vol += d_JxW[iG];
  
  return vol;
  
}


