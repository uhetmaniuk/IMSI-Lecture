#include "Elements/fe1DP1.h"
#include <iostream>

using namespace std;
using namespace FECode;


fe1DP1::fe1DP1
(
 const std::vector<double> &nodeCoord,
 const std::vector<double> &lQ,
 const std::vector<double> &wQ
)
: 
d_JxW(vector<double>(wQ.size())),
d_numGauss(wQ.size()),
d_Phi(std::vector<double>(2*wQ.size())),
d_PhiGrad(std::vector<double>(2*wQ.size())),
d_PointInt(std::vector<double>(wQ.size())),
d_x(std::vector<double>(2)),
d_lQuad(lQ),
d_wQuad(wQ)
{

  d_x[0] = nodeCoord[0]; d_x[1] = nodeCoord[1];

  double h = d_x[1] - d_x[0];
  double xM = (d_x[1] + d_x[0])*0.5;

  double xi;
  for (size_t ii = 0; ii < d_numGauss; ++ii)
  {
    xi = d_lQuad[ii];
    //--------
    d_Phi[ii] = (1.0 - xi) * 0.5;
    d_PhiGrad[ii] = -1.0 / h;
    //--------
    d_Phi[ii + d_numGauss] = (1.0 + xi) * 0.5;
    d_PhiGrad[ii + d_numGauss] = 1.0 / h;
    //--------
    d_JxW[ii] = d_wQuad[ii] * 0.5 * h;
    d_PointInt[ii] = xM + xi * h;
  }
  
}


void fe1DP1::GetShapeGrad
(
  const size_t iP
  , const size_t iQ
  , std::valarray<double> &sg
)
{

  if (sg.size() != 1)
    sg.resize(1);
  sg[0] = d_PhiGrad[iP * d_numGauss + iQ];

}


void fe1DP1::GetShapeHessian
(
  const size_t iP
  , const size_t iQ
  , std::vector<double> &sh
)
{

  sh = std::vector<double>(1, 0.0);

}


double fe1DP1::Volume
() 
{

  return d_x[1] - d_x[0];

}


/*
vector<double> fe1DP1::NodalGrad(unsigned int iP) {
//
// Extrapolate the gradient from the value at Gauss points
//

  vector<double> res(1, 0.0);

  double coeff = 0.5*(1.0 + sqrt(3.0));

  double dphiR = coeff*PhiGrad[2*iP+1] + (1.0-coeff)*PhiGrad[2*iP];
  double dphiL = (1.0-coeff)*PhiGrad[2*iP+1] + coeff*PhiGrad[2*iP];

  res[0] = x[0]*dphiR + x[1]*dphiL;

  return res;

}
*/


