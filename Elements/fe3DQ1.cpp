#include "Elements/fe3DQ1.h"

#include <cassert>
#include <cmath>
#include <iostream>


using namespace std;
using namespace FECode;


fe3DQ1::fe3DQ1
(
  const std::vector<double> &nodeCoord,
  const std::vector<double> &lQ,
  const std::vector<double> &wQ
)
: 
  d_JxW(vector<double>(wQ.size())),
  d_numGauss(wQ.size()),
  d_Phi(std::vector<double>(nnode*wQ.size())),
  d_PhiGrad(std::vector<double>(sdim*nnode*wQ.size())),
  d_PointInt(std::vector<double>(sdim*wQ.size())),
  d_lQuad(lQ),
  d_wQuad(wQ)
{

  std::vector<double> x(nnode), y(nnode), z(nnode);
  for (int j = 0; j < nnode; ++j)
  {
    x[j] = nodeCoord[sdim*j];
    y[j] = nodeCoord[sdim*j+1];
    z[j] = nodeCoord[sdim*j+2];
  }

  double xi, eta, zeta;
  double J[9], Jac = 0.0;
  double dqxi[8], dqeta[8], dqzeta[8];

  for (size_t ipp = 0; ipp < d_numGauss; ++ipp)
  {
    xi = d_lQuad[sdim*ipp];
    eta = d_lQuad[sdim*ipp + 1];
    zeta = d_lQuad[sdim*ipp + 2];
    //---
    double halfOnePXi = 0.5*(1.0 + xi);
    double halfOneMXi = 1.0 - halfOnePXi;
    double halfOnePEta = 0.5*(1.0 + eta);
    double halfOneMEta = 1.0 - halfOnePEta;
    double halfOnePZeta = 0.5*(1.0 + zeta);
    double halfOneMZeta = 1.0 - halfOnePZeta;
    //---------------------
    size_t pos = ipp;
    d_Phi[pos] = halfOneMXi * halfOneMEta * halfOneMZeta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOneMEta * halfOneMZeta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOnePEta * halfOneMZeta; pos += d_numGauss;
    d_Phi[pos] = halfOneMXi * halfOnePEta * halfOneMZeta; pos += d_numGauss;
    d_Phi[pos] = halfOneMXi * halfOneMEta * halfOnePZeta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOneMEta * halfOnePZeta; pos += d_numGauss;
    d_Phi[pos] = halfOnePXi * halfOnePEta * halfOnePZeta; pos += d_numGauss;
    d_Phi[pos] = halfOneMXi * halfOnePEta * halfOnePZeta;
    //---------------------
    dqxi[0] = -0.5 * halfOneMEta * halfOneMZeta;
    dqxi[1] = 0.5 * halfOneMEta * halfOneMZeta;
    dqxi[2] = 0.5 * halfOnePEta * halfOneMZeta;
    dqxi[3] = -0.5 * halfOnePEta * halfOneMZeta;
    dqxi[4] = -0.5 * halfOneMEta * halfOnePZeta;
    dqxi[5] = 0.5 * halfOneMEta * halfOnePZeta;
    dqxi[6] = 0.5 * halfOnePEta * halfOnePZeta;
    dqxi[7] = -0.5 * halfOnePEta * halfOnePZeta;
    //---
    dqeta[0] = halfOneMXi * (-0.5) * halfOneMZeta;
    dqeta[1] = halfOnePXi * (-0.5) * halfOneMZeta;
    dqeta[2] = halfOnePXi * (0.5) * halfOneMZeta;
    dqeta[3] = halfOneMXi * (0.5) * halfOneMZeta;
    dqeta[4] = halfOneMXi * (-0.5) * halfOnePZeta;
    dqeta[5] = halfOnePXi * (-0.5) * halfOnePZeta;
    dqeta[6] = halfOnePXi * (0.5) * halfOnePZeta;
    dqeta[7] = halfOneMXi * (0.5) * halfOnePZeta;
    //---
    dqzeta[0] = halfOneMXi * halfOneMEta * (-0.5);
    dqzeta[1] = halfOnePXi * halfOneMEta * (-0.5);
    dqzeta[2] = halfOnePXi * halfOnePEta * (-0.5);
    dqzeta[3] = halfOneMXi * halfOnePEta * (-0.5);
    dqzeta[4] = halfOneMXi * halfOneMEta * (0.5);
    dqzeta[5] = halfOnePXi * halfOneMEta * (0.5);
    dqzeta[6] = halfOnePXi * halfOnePEta * (0.5);
    dqzeta[7] = halfOneMXi * halfOnePEta * (0.5);
    //---
    J[0] = 0.0; J[3] = 0.0; J[6] = 0.0;
    J[1] = 0.0; J[4] = 0.0; J[7] = 0.0;
    J[2] = 0.0; J[5] = 0.0; J[8] = 0.0;
    for (size_t j = 0; j < nnode; ++j)
    {
      J[0] += x[j]*dqxi[j]; J[3] += y[j]*dqxi[j]; J[6] += z[j]*dqxi[j];
      J[1] += x[j]*dqeta[j]; J[4] += y[j]*dqeta[j]; J[7] += z[j]*dqeta[j];
      J[2] += x[j]*dqzeta[j]; J[5] += y[j]*dqzeta[j]; J[8] += z[j]*dqzeta[j];
    }
    //---------------------
    inverseJacobian(sdim, &(J[0]), Jac);
    if (Jac <= 0.0) 
    {
      if (Jac == 0.0)
        cerr << "\n fe3DQ1 >> Zero Jacobian determinant !!!\n\n";
      else
        cerr << "\n fe3DQ1 >> Negative Jacobian determinant !!!\n\n";
      bool hasPositiveJacobian = false;
      assert(hasPositiveJacobian == true);
    }
    //---------------------
    d_JxW[ipp] = Jac * d_wQuad[ipp];
    //---------------------
    for (unsigned int j = 0; j < nnode; ++j)
    {
      pos = sdim * d_numGauss * j + sdim * ipp;
      d_PhiGrad[pos]   = J[0]*dqxi[j] + J[3]*dqeta[j] + J[6]*dqzeta[j];
      d_PhiGrad[pos+1] = J[1]*dqxi[j] + J[4]*dqeta[j] + J[7]*dqzeta[j];
      d_PhiGrad[pos+2] = J[2]*dqxi[j] + J[5]*dqeta[j] + J[8]*dqzeta[j];
    }
    //---------------------
    double tmpX = 0.0, tmpY = 0.0, tmpZ = 0.0;
    for (size_t j = 0; j < nnode; ++j)
    {
      double shape = d_Phi[j * d_numGauss + ipp];
      tmpX += x[j]*shape;
      tmpY += y[j]*shape;
      tmpZ += z[j]*shape;
    }
    pos = sdim * ipp;
    d_PointInt[pos] = tmpX;
    d_PointInt[pos + 1] = tmpY;
    d_PointInt[pos + 2] = tmpZ;
    //---------------------

  } // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

}


void fe3DQ1::GetShapeGrad
(
  const size_t iP
  , const size_t iQ
  , std::valarray<double> &sg
)
{

  size_t pos = sdim * d_numGauss * iP + sdim * iQ;
  if (sg.size() != sdim)
    sg.resize(sdim);
  sg[0] = d_PhiGrad[pos];
  sg[1] = d_PhiGrad[pos + 1];
  sg[2] = d_PhiGrad[pos + 2];

}


void fe3DQ1::GetShapeHessian
(
  const size_t iP
  , const size_t iQ
  , std::vector<double> &sh
) 
{

  sh = std::vector<double>(8);

}


double fe3DQ1::Volume
() 
{

  double vol = 0.0;
  for (int iG = 0; iG < d_numGauss; ++iG)
  {
    vol += d_JxW[iG];
  }

  return vol;

}


