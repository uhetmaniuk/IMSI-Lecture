#include "Elements/LagrangeElement.h"

#include <cassert>
#include <iostream>


void FECode::LagrangeElement::inverseJacobian
(
 int n
 , double *J
 , double &Jac
) 
{

// This routine computes the inverse of a matrix n x n stored in J.
// The input matrix should be stored column-wise.
// The inverse is then stored in J.
// Jac returns the determinant of J.

  if (n == 1) 
  {
    Jac = J[0];
    J[0] = 1.0/J[0];
  }
  else if (n == 2) 
  {
    double J11 = J[0];
    Jac = J11*J[3] - J[1]*J[2];
    double invJac = 1.0/Jac;
    J[0] = J[3] * invJac;
    J[1] = -J[1] * invJac;
    J[2] = -J[2] * invJac;
    J[3] = J11 * invJac;
  }
  else if (n == 3) 
  {
    double iJ11, iJ12, iJ13, iJ21, iJ22, iJ23, iJ31, iJ32, iJ33;
    iJ11 = J[4]*J[8] - J[5]*J[7];
    iJ12 = J[6]*J[5] - J[8]*J[3];
    iJ13 = J[3]*J[7] - J[4]*J[6];
    iJ21 = J[7]*J[2] - J[1]*J[8];
    iJ22 = J[0]*J[8] - J[2]*J[6];
    iJ23 = J[6]*J[1] - J[7]*J[0];
    iJ31 = J[1]*J[5] - J[2]*J[4];
    iJ32 = J[3]*J[2] - J[5]*J[0];
    iJ33 = J[0]*J[4] - J[1]*J[3];
    //-------------
    Jac = J[0]*iJ11 + J[3]*iJ21 + J[6]*iJ31;
    //-------------
    double invJac = 1.0/Jac;
    J[0] = iJ11 * invJac; J[3] = iJ12 * invJac; J[6] = iJ13 * invJac;
    J[1] = iJ21 * invJac; J[4] = iJ22 * invJac; J[7] = iJ23 * invJac;
    J[2] = iJ31 * invJac; J[5] = iJ32 * invJac; J[8] = iJ33 * invJac;
  }
  else {
    // TO BE DONE
  }

  if (Jac == 0.0) 
  {
    std::cerr << "\n LagrangeElement >> Zero Jacobian determinant !!!\n\n";
    bool hasPositiveJacobian = false;
    assert(hasPositiveJacobian == true);
  }
  
}

