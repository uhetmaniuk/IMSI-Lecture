#include "Elements/fe2DP1.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>


using namespace std;
using namespace FECode;


fe2DP1::fe2DP1
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

}


void fe2DP1::makeShapeValue
		(
		)
{

	double xi, eta;

	double y12 = d_y[1] - d_y[2];
	double x12 = d_x[1] - d_x[2];

	double det = d_x[0] * y12 - d_y[0] * x12 + d_x[1] * d_y[2] - d_x[2] * d_y[1];
	if (det <= DBL_MIN) {
		if (det == 0.0)
			cerr << "\n fe2DP1 >> Zero Jacobian determinant !!!\n\n";
		else
			cerr << "\n fe2DP1 >> Negative Jacobian determinant !!!\n\n";
		bool hasPositiveJacobian = false;
		assert(hasPositiveJacobian == true);
	}

	double cdet = 1.0 / det;

	d_PhiGrad[0] = cdet * y12;
	d_PhiGrad[1] = cdet * (-x12);

	d_PhiGrad[sdim*d_numGauss] = cdet * (d_y[2] - d_y[0]);
	d_PhiGrad[sdim*d_numGauss + 1] = cdet * (d_x[0] - d_x[2]);

	d_PhiGrad[2*sdim*d_numGauss] = cdet * (d_y[0] - d_y[1]);
	d_PhiGrad[2*sdim*d_numGauss + 1] = cdet * (d_x[1] - d_x[0]);

	for (size_t ipp = 0; ipp < d_numGauss; ++ipp) {
		xi = d_lQuad[sdim*ipp];
		eta = d_lQuad[sdim*ipp + 1];
		//--------
		unsigned int pos = ipp;
		d_Phi[pos] = xi; pos += d_numGauss;
		d_Phi[pos] = eta; pos += d_numGauss;
		d_Phi[pos] = 1.0 - xi - eta;
		//---------------------
		d_JxW[ipp] = d_wQuad[ipp] * det * 0.5;
		//---------------------
		for (unsigned int j = 0; j < nnode; ++j) {
			pos = sdim*d_numGauss*j + sdim*ipp;
			d_PhiGrad[pos]   = d_PhiGrad[sdim*d_numGauss*j];
			d_PhiGrad[pos+1] = d_PhiGrad[sdim*d_numGauss*j + 1];
		}
		//---------------------

	} // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

	d_hasShapeValue = true;

}


void fe2DP1::makePointInt
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



void fe2DP1::GetShapeGrad
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


void fe2DP1::GetShapeHessian
		(
				const size_t iP
				, const size_t iQ
				, std::vector<double> &sh
		)
{

	for (auto &val : sh)
		val = 0.0;

}


double fe2DP1::Volume
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


