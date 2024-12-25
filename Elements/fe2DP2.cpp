#include "Elements/fe2DP2.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>


using namespace std;
using namespace FECode;


fe2DP2::fe2DP2
		(
				const std::vector<double> &nodeCoord,
				const std::vector<double> &lQ,
				const std::vector<double> &wQ
		)
		:
		d_JxW(wQ.size()),
		d_numGauss(wQ.size()),
		d_Phi(nnode*wQ.size()),
		d_PhiGrad(sdim*nnode*wQ.size()),
		d_PointInt(sdim*wQ.size()),
		d_x(nnode),
		d_y(nnode),
		d_lQuad(lQ),
		d_wQuad(wQ),
		d_hasShapeValue(false),
		d_hasPointInt(false)
{

	d_x[0] = nodeCoord[0]; d_y[0] = nodeCoord[1];
	d_x[1] = nodeCoord[2]; d_y[1] = nodeCoord[3];
	d_x[2] = nodeCoord[4]; d_y[2] = nodeCoord[5];
	d_x[3] = nodeCoord[6]; d_y[3] = nodeCoord[7];
	d_x[4] = nodeCoord[8]; d_y[4] = nodeCoord[9];
	d_x[5] = nodeCoord[10]; d_y[5] = nodeCoord[11];

}


void fe2DP2::makeShapeValue
		(
		)
{
//
// This code is translated from TRIG6SHAPE
// 

	double zeta0, zeta1, zeta2;
	double jx0, jx1, jx2, jy0, jy1, jy2;
	double det, cdet;
	double x10, x21, x02, y01, y12, y20;

	for (size_t ipp = 0; ipp < d_numGauss; ++ipp) {

		zeta0 = d_lQuad[2*ipp];
		zeta1 = d_lQuad[2*ipp + 1];
		zeta2 = 1.0 - zeta0 - zeta1;

		d_Phi[ipp]                = zeta0 * (2.0 * zeta0 - 1.0);
		d_Phi[ipp +   d_numGauss] = zeta1 * (2.0 * zeta1 - 1.0);
		d_Phi[ipp + 2*d_numGauss] = zeta2 * (2.0 * zeta2 - 1.0);
		d_Phi[ipp + 3*d_numGauss] = 4.0 * zeta0 * zeta1;
		d_Phi[ipp + 4*d_numGauss] = 4.0 * zeta1 * zeta2;
		d_Phi[ipp + 5*d_numGauss] = 4.0 * zeta2 * zeta0;

		jx0 = 4.0 * (d_x[0]*zeta0 + d_x[3]*zeta1 + d_x[5]*zeta2) - d_x[0];
		jx1 = 4.0 * (d_x[1]*zeta1 + d_x[4]*zeta2 + d_x[3]*zeta0) - d_x[1];
		jx2 = 4.0 * (d_x[2]*zeta2 + d_x[5]*zeta0 + d_x[4]*zeta1) - d_x[2];

		jy0 = 4.0 * (d_y[0]*zeta0 + d_y[3]*zeta1 + d_y[5]*zeta2) - d_y[0];
		jy1 = 4.0 * (d_y[1]*zeta1 + d_y[4]*zeta2 + d_y[3]*zeta0) - d_y[1];
		jy2 = 4.0 * (d_y[2]*zeta2 + d_y[5]*zeta0 + d_y[4]*zeta1) - d_y[2];

		x10 = jx1 - jx0;
		x21 = jx2 - jx1;
		x02 = jx0 - jx2;

		y01 = jy0 - jy1;
		y12 = jy1 - jy2;
		y20 = jy2 - jy0;

		det = x10 * y20 - y01 * x02;
		if (det <= DBL_MIN) {
	if (det == 0.0)
		cerr << "\n fe2DP2 >> Zero Jacobian determinant !!!\n\n";
	else
		cerr << "\n fe2DP2 >> Negative Jacobian determinant !!!\n\n";
	bool hasPositiveJacobian = false;
	assert(hasPositiveJacobian == true);
		}

		d_JxW[ipp] = d_wQuad[ipp] * det * 0.5;

		cdet = 1.0 / det;

		//--- x-component of gradients for nodes 0, 1, 2
		 d_PhiGrad[ipp*sdim]                     = cdet * (4.0*zeta0 - 1.0) * y12;
		 d_PhiGrad[ipp*sdim + sdim*d_numGauss]   = cdet * (4.0*zeta1 - 1.0) * y20;
		 d_PhiGrad[ipp*sdim + sdim*d_numGauss*2] = cdet * (4.0*zeta2 - 1.0) * y01;

		//--- y-component of gradients for nodes 0, 1, 2
		 d_PhiGrad[ipp*sdim + 1]                     = cdet * (4.0*zeta0 - 1.0) * x21;
		 d_PhiGrad[ipp*sdim + 1 + sdim*d_numGauss]   = cdet * (4.0*zeta1 - 1.0) * x02;
		 d_PhiGrad[ipp*sdim + 1 + sdim*d_numGauss*2] = cdet * (4.0*zeta2 - 1.0) * x10;

		//--- x-component of gradients for nodes 3, 4, 5
		 d_PhiGrad[ipp*sdim + sdim*d_numGauss*3] = cdet*4.0*(zeta1*y12 + zeta0*y20);
		 d_PhiGrad[ipp*sdim + sdim*d_numGauss*4] = cdet*4.0*(zeta2*y20 + zeta1*y01);
		 d_PhiGrad[ipp*sdim + sdim*d_numGauss*5] = cdet*4.0*(zeta0*y01 + zeta2*y12);

		//--- y-component of gradients for nodes 3, 4, 5
		 d_PhiGrad[ipp*sdim + 1 + sdim*d_numGauss*3] = cdet*4.0*(zeta1*x21 + zeta0*x02);
		 d_PhiGrad[ipp*sdim + 1 + sdim*d_numGauss*4] = cdet*4.0*(zeta2*x02 + zeta1*x10);
		 d_PhiGrad[ipp*sdim + 1 + sdim*d_numGauss*5] = cdet*4.0*(zeta0*x10 + zeta2*x21);

	} // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

	d_hasShapeValue = true;

}


void fe2DP2::makePointInt
		(
		)
{

	if (d_hasShapeValue == false)
		makeShapeValue();

	unsigned int iN, iQ;
	for (iQ = 0; iQ < d_numGauss; ++iQ) {
		double tmpX = 0.0, tmpY = 0.0;
		for (iN = 0; iN < nnode; ++iN) {
			double shape = d_Phi[iN*d_numGauss + iQ];
			tmpX += d_x[iN]*shape;
			tmpY += d_y[iN]*shape;
		}
		unsigned int pos = sdim*iQ;
		d_PointInt[pos] = tmpX;
		d_PointInt[pos+1] = tmpY;
	}

	d_hasPointInt = true;

}



void fe2DP2::GetShapeGrad
		(
				const size_t iP
				, const size_t iQ
				, std::valarray<double> &sg
		)
{

	if (d_hasShapeValue == false)
		makeShapeValue();

	size_t pos = sdim * d_numGauss * iP + sdim * iQ;
	if (sg.size() != sdim)
		sg.resize(sdim);
	sg[0] = d_PhiGrad[pos];
	sg[1] = d_PhiGrad[pos+1];

}


void fe2DP2::GetShapeHessian
		(
				const size_t iP
				, const size_t iQ
				, std::vector<double> &sh
		)
{

	/* --- 2019/08 = Not Implemented */

	sh[0] = 0.0;
	sh[1] = 0.0;
	sh[2] = sh[1];
	sh[3] = 0.0;

}


double fe2DP2::Volume
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


