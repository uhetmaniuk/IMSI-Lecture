#include "Elements/fe3DP2.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>


using namespace std;
using namespace FECode;


fe3DP2::fe3DP2
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

	for (int j = 0; j < nnode; ++j) {
		d_x[j] = nodeCoord[j*sdim];
		d_y[j] = nodeCoord[j*sdim + 1];
		d_z[j] = nodeCoord[j*sdim + 2];
	}

}


void fe3DP2::makeShapeValue
		(
		)
{

	double J[9];
	double Jac = 0.0;
	double zeta0, zeta1, zeta2, zeta3;

	double dq0[10], dq1[10], dq2[10];
	for (int ii = 0; ii < nnode; ++ii) {
		dq0[ii] = 0.0;
		dq1[ii] = 0.0;
		dq2[ii] = 0.0;
	}

	for (size_t ipp = 0; ipp < d_numGauss; ++ipp) {

		zeta0 = d_lQuad[sdim*ipp];
		zeta1 = d_lQuad[sdim*ipp + 1];
		zeta2 = d_lQuad[sdim*ipp + 2];
		zeta3 = 1.0 - zeta0 - zeta1 - zeta2;

		d_Phi[ipp]                = zeta0*(2.0*zeta0 - 1.0);
		d_Phi[ipp + d_numGauss]   = zeta1*(2.0*zeta1 - 1.0);
		d_Phi[ipp + 2*d_numGauss] = zeta2*(2.0*zeta2 - 1.0);
		d_Phi[ipp + 3*d_numGauss] = zeta3*(2.0*zeta3 - 1.0);

		d_Phi[ipp + 4*d_numGauss] = 4.0*zeta0*zeta1;
		d_Phi[ipp + 5*d_numGauss] = 4.0*zeta0*zeta2;
		d_Phi[ipp + 6*d_numGauss] = 4.0*zeta0*zeta3;
		d_Phi[ipp + 7*d_numGauss] = 4.0*zeta1*zeta2;
		d_Phi[ipp + 8*d_numGauss] = 4.0*zeta2*zeta3;
		d_Phi[ipp + 9*d_numGauss] = 4.0*zeta3*zeta1;

		dq0[0] = 4*zeta0 - 1.0;
		dq0[3] = 1.0 - 4*zeta3;

		dq0[4] = 4*zeta1;
		dq0[5] = 4*zeta2;
		dq0[6] = 4*zeta3 - 4.0*zeta0;

		dq0[8] = -4.0*zeta2;
		dq0[9] = -4.0*zeta1;

		dq1[1] = 4*zeta1 - 1.0;
		dq1[3] = 1.0 - 4*zeta3;

		dq1[4] = 4*zeta0;
		dq1[6] = -4*zeta0;
		dq1[7] = 4*zeta2;

		dq1[8] = -4.0*zeta2;
		dq1[9] = 4.0*zeta3 - 4.0*zeta1;

		dq2[2] = 4*zeta2 - 1.0;
		dq2[3] = 1.0 - 4*zeta3;

		dq2[5] = 4*zeta0;
		dq2[6] = -4*zeta0;
		dq2[7] = 4*zeta1;

		dq2[8] = 4.0*zeta3 - 4.0*zeta2;
		dq2[9] = -4.0*zeta1;

		for (int jj = 0; jj < 9; ++jj)
			J[jj] = 0.0;

		for (int j = 0; j < nnode; ++j){
			J[0] += d_x[j] * dq0[j];
			J[3] += d_y[j] * dq0[j];
			J[6] += d_z[j] * dq0[j];

			J[1] += d_x[j] * dq1[j];
			J[4] += d_y[j] * dq1[j];
			J[7] += d_z[j] * dq1[j];

			J[2] += d_x[j] * dq2[j];
			J[5] += d_y[j] * dq2[j];
			J[8] += d_z[j] * dq2[j];
		}

		inverseJacobian(3, &(J[0]), Jac);

		if (Jac < DBL_MIN) {
			if (Jac == 0.0)
				cerr << "\n fe3DP2 >> Zero Jacobian determinant !!!\n\n";
			else
				cerr << "\n fe3DP2 >> Negative Jacobian determinant !!!\n\n";
			bool hasPositiveJacobian = false;
			assert(hasPositiveJacobian == true);
		}

		d_JxW[ipp] = d_wQuad[ipp] * Jac * 0.5;

		for (int j = 0; j < nnode; ++j)
		{
			int pos = sdim * d_numGauss * j + sdim * ipp;
			d_PhiGrad[pos]   = J[0]*dq0[j] + J[3]*dq1[j] + J[6]*dq2[j];
			d_PhiGrad[pos+1] = J[1]*dq0[j] + J[4]*dq1[j] + J[7]*dq2[j];
			d_PhiGrad[pos+2] = J[2]*dq0[j] + J[5]*dq1[j] + J[8]*dq2[j];
		}

	} // for (size_t ipp = 0; ipp < d_numGauss; ++ipp)

	d_hasShapeValue = true;

}


void fe3DP2::makePointInt
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



void fe3DP2::GetShapeGrad
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


void fe3DP2::GetShapeHessian
		(
				const size_t iP
				, const size_t iQ
				, std::vector<double> &sh
		)
{

	/* --- Not Implemented Yet --- */
	assert(0 > 1);
	exit(-1);

}


double fe3DP2::Volume
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


