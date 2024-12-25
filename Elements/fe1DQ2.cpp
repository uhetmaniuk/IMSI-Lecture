#include "Elements/fe1DQ2.h"


using namespace std;
using namespace FECode;


fe1DQ2::fe1DQ2
(
 const std::vector<double> &nodeCoord,
 const std::vector<double> &lQ,
 const std::vector<double> &wQ
)
: 
d_JxW(vector<double>(wQ.size())),
d_numGauss(wQ.size()),
d_Phi(std::vector<double>(3*wQ.size())),
d_PhiGrad(std::vector<double>(3*wQ.size())),
d_PointInt(std::vector<double>(wQ.size())),
d_x(std::vector<double>(3)),
d_lQuad(lQ),
d_wQuad(wQ)
{

  d_x[0] = nodeCoord[0];
  d_x[1] = nodeCoord[1];
  d_x[2] = nodeCoord[2];

  double xi;
  double J;

  for (size_t iQ = 0; iQ < d_numGauss; ++iQ) {
    xi = d_lQuad[iQ];
    d_Phi[iQ] = 0.5*(xi - 1.0)*xi;
    d_Phi[d_numGauss + iQ] = 0.5*(1.0 + xi)*xi;
    d_Phi[2*d_numGauss + iQ] = 1.0 - xi*xi;
    d_PhiGrad[iQ] = xi - 0.5;
    d_PhiGrad[d_numGauss + iQ] = xi + 0.5;
    d_PhiGrad[2*d_numGauss + iQ] = -2.0*xi;
    J = d_PhiGrad[iQ]*d_x[0] + d_PhiGrad[d_numGauss+iQ]*d_x[1] + d_PhiGrad[2*d_numGauss+iQ]*d_x[2];
    d_PhiGrad[iQ] = d_PhiGrad[iQ]/J;
    d_PhiGrad[d_numGauss + iQ] = d_PhiGrad[d_numGauss + iQ]/J;
    d_PhiGrad[2*d_numGauss + iQ] = d_PhiGrad[2*d_numGauss + iQ]/J;
    d_JxW[iQ] = d_wQuad[iQ]*J;
  }

  for (size_t i = 0; i < d_numGauss; ++i)
  {
    d_PointInt[i] = d_x[0]*d_Phi[i] + d_x[1]*d_Phi[d_numGauss+i]
                    + d_x[2]*d_Phi[2*d_numGauss+i];
  }

}


void fe1DQ2::GetShapeGrad
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


void fe1DQ2::GetShapeHessian
(
  const size_t iP
  , const size_t iQ
  , std::vector<double> &sh
)
{

  sh = std::vector<double>(1, 0.0);

  double d2phidxi2[3];
  d2phidxi2[0] = 1.0;
  d2phidxi2[1] = 1.0;
  d2phidxi2[2] = -2.0;

  double d2xdxi2 = d_x[0]*d2phidxi2[0] + d_x[1]*d2phidxi2[1] + d_x[2]*d2phidxi2[2];
  double oneOverJ = d_wQuad[iQ]/d_JxW[iQ];

  sh[0] = d2phidxi2[iP] - d_PhiGrad[iP*d_numGauss + iQ]*d2xdxi2;
  sh[0] = sh[0]*(oneOverJ*oneOverJ);

}


double fe1DQ2::Volume() 
{

  return d_x[1] - d_x[0];

} 
  

/*
vector<double> fe1DQ2::NodalGrad(unsigned int iP) {
//
// Extrapolate the gradient from the value at Gauss points
//

  vector<double> res(1, 0.0);

  double dphi[3];
  
  double coeff[3];
  
  unsigned int pos = numGauss*iP;
  dphi[i] = coeff[0]*PhiGrad[pos] + coeff[1]*PhiGrad[pos+1] + coeff[2]*PhiGrad[pos+2];
  double dphiL = (1.0-coeff)*PhiGrad[pos+1] + coeff*PhiGrad[pos] + 0.0*PhiGrad[pos+2];
  double dphiC = 0.0;

  res[0] = x[0]*dphiR + x[1]*dphiL + x[2]*dphiC;

  return res;

}
*/


