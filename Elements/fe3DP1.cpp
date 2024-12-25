#include "Elements/fe3DP1.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>


using namespace std;
using namespace FECode;


fe3DP1::fe3DP1
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
		d_z(nnode),
		d_lQuad(lQ),
		d_wQuad(wQ),
		d_hasShapeValue(false),
		d_hasPointInt(false)
{

	d_x[0] = nodeCoord[0]; d_y[0] = nodeCoord[1]; d_z[0] = nodeCoord[2];
	d_x[1] = nodeCoord[3]; d_y[1] = nodeCoord[4]; d_z[1] = nodeCoord[5];
	d_x[2] = nodeCoord[6]; d_y[2] = nodeCoord[7]; d_z[2] = nodeCoord[8];
	d_x[3] = nodeCoord[9]; d_y[3] = nodeCoord[10]; d_z[3] = nodeCoord[11];

}


void fe3DP1::makeShapeValue
		(
		)
{

	double J[9];

	J[0] = d_x[0] - d_x[3];
	J[1] = d_x[1] - d_x[3];
	J[2] = d_x[2] - d_x[3];

	J[3] = d_y[0] - d_y[3];
	J[4] = d_y[1] - d_y[3];
	J[5] = d_y[2] - d_y[3];

	J[6] = d_z[0] - d_z[3];
	J[7] = d_z[1] - d_z[3];
	J[8] = d_z[2] - d_z[3];

	double Jac = 0.0;
	inverseJacobian(3, &(J[0]), Jac);

	if (Jac < DBL_MIN) {
		if (Jac == 0.0)
			cerr << "\n fe3DP1 >> Zero Jacobian determinant !!!\n\n";
		else
			cerr << "\n fe3DP1 >> Negative Jacobian determinant !!!\n\n";
		bool hasPositiveJacobian = false;
		assert(hasPositiveJacobian == true);
	}

	// --- Derivative is [1; 0; 0]
	d_PhiGrad[0] = J[0];
	d_PhiGrad[1] = J[1];
	d_PhiGrad[2] = J[2];

	// --- Derivative is [0; 1; 0]
	d_PhiGrad[sdim*d_numGauss + 0] = J[3];
	d_PhiGrad[sdim*d_numGauss + 1] = J[4];
	d_PhiGrad[sdim*d_numGauss + 2] = J[5];

	// --- Derivative is [0; 0; 1]
	d_PhiGrad[sdim*d_numGauss*2 + 0] = J[6];
	d_PhiGrad[sdim*d_numGauss*2 + 1] = J[7];
	d_PhiGrad[sdim*d_numGauss*2 + 2] = J[8];

	// --- Derivative is [-1; -1; -1]
	d_PhiGrad[sdim*d_numGauss*3 + 0] = -(J[0] + J[3] + J[6]);
	d_PhiGrad[sdim*d_numGauss*3 + 1] = -(J[1] + J[4] + J[7]);
	d_PhiGrad[sdim*d_numGauss*3 + 2] = -(J[2] + J[5] + J[8]);

	double zeta0, zeta1, zeta2;

	for (size_t ipp = 0; ipp < d_numGauss; ++ipp) {

		zeta0 = d_lQuad[sdim*ipp];
		zeta1 = d_lQuad[sdim*ipp + 1];
		zeta2 = d_lQuad[sdim*ipp + 2];

		d_Phi[ipp]                = zeta0;
		d_Phi[ipp + d_numGauss]   = zeta1;
		d_Phi[ipp + 2*d_numGauss] = zeta2;
		d_Phi[ipp + 3*d_numGauss] = 1.0 - zeta0 - zeta1 - zeta2;

		d_JxW[ipp] = d_wQuad[ipp] * Jac * 0.5;

		for (unsigned int j = 0; j < nnode; ++j) {
			unsigned int pos = sdim*d_numGauss*j + sdim*ipp;
			d_PhiGrad[pos]   = d_PhiGrad[sdim*d_numGauss*j];
			d_PhiGrad[pos+1] = d_PhiGrad[sdim*d_numGauss*j + 1];
			d_PhiGrad[pos+2] = d_PhiGrad[sdim*d_numGauss*j + 2];
		}

	} // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

	d_hasShapeValue = true;

}


void fe3DP1::makePointInt
		(
		)
{

	if (!d_hasShapeValue)
		makeShapeValue();

	unsigned int iN, iQ;
	for (iQ = 0; iQ < d_numGauss; ++iQ) {
		double tmpX = 0.0, tmpY = 0.0, tmpZ = 0.0;
		for (iN = 0; iN < nnode; ++iN) {
			double shape = d_Phi[iN*d_numGauss + iQ];
			tmpX += d_x[iN]*shape;
			tmpY += d_y[iN]*shape;
			tmpZ += d_z[iN]*shape;
		}
		unsigned int pos = sdim*iQ;
		d_PointInt[pos] = tmpX;
		d_PointInt[pos+1] = tmpY;
		d_PointInt[pos+2] = tmpZ;
	}

	d_hasPointInt = true;

}



void fe3DP1::GetShapeGrad
		(
				const size_t iP
				, const size_t iQ
				, std::valarray<double> &sg
		)
{

	if (!d_hasShapeValue)
		makeShapeValue();

	size_t pos = sdim * d_numGauss * iP + sdim * iQ;
	if (sg.size() != sdim)
		sg.resize(sdim);
	sg[0] = d_PhiGrad[pos];
	sg[1] = d_PhiGrad[pos+1];
	sg[2] = d_PhiGrad[pos+2];

}


void fe3DP1::GetShapeHessian
		(
				const size_t iP
				, const size_t iQ
				, std::vector<double> &sh
		)
{

	for (auto &val : sh)
		val = 0.0;

}


double fe3DP1::Volume
		(
		)
{

	if (!d_hasShapeValue)
		makeShapeValue();

	double vol = 0.0;
	for (int iG = 0; iG < d_numGauss; ++iG)
		vol += d_JxW[iG];

	return vol;

}


