#include "Geometry/Hexahedron.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "Geometry/Quadrilateral.h"
#include "Geometry/Vertex.h"
#include "Math/MathCXXCore.hpp"


using namespace FECode;


Hexahedron::Hexahedron
(
  short dim, 
  std::vector<Vertex*> &vertexList,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, HEXAHEDRON)
{

  // Test whether the dimension is at least 3
  if (dim < 3) {
    std::cerr << "\n !!! An hexahedron requires three dimensions !!!\n\n";
    assert(dim >= 3);
  }

  // Test whether the vertex list is of size 8.
  if (vertexList.size() != 8) {
    std::cerr << "\n !!! An hexahedron requires eight distinct vertices !!!\n\n";
    assert(vertexList.size() == 8);
  }

  // Copy the vertex list.
  d_vertexList.assign(vertexList.begin(), vertexList.end());

}


Hexahedron::~Hexahedron() {
}


std::vector<double> Hexahedron::GetCoordinatesList() const
{

  std::vector<double> coord(d_vertexList.size() * d_dim);

  int count = 0;
  for (int ii = 0; ii < d_vertexList.size(); ++ii)
  {
    std::vector<double> node = d_vertexList[ii]->GetCoordinates();
    for (int jj = 0; jj < node.size(); ++jj)
    {
      coord[count] = node[jj];
      count += 1;
    }
  }

  return coord;

}


bool Hexahedron::GetMapJacobian
(
  const std::vector<double> &lCoord,
  std::vector<double> &B
) const 
{

  if (B.size() != 3*d_dim)
  {
    B.resize(3*d_dim);
  }

  if (d_dim > 3) {
    std::cerr << "\n !! Hexahedron::GetMapJacobian >> ";
    std::cerr << "The Jacobian in high dimensions is not implemented !!\n\n";
    assert(0 > 1);
    return false;
  }  

  if (d_dim == 3) {
    //
    // Compute the Jacobian for an hexahedron
    //
    double xi = lCoord[0];
    double eta = lCoord[1];
    double zeta = lCoord[2];
    //---------------------
    std::vector<double> dqxi(8);
    dqxi[1] = (1.0 - eta) * (1.0 - zeta);
    dqxi[2] = eta * (1.0 - zeta);
    dqxi[5] = (1.0 - eta) * zeta; 
    dqxi[6] = eta * zeta; 
    dqxi[0] = -dqxi[1];
    dqxi[3] = -dqxi[2];
    dqxi[4] = -dqxi[5]; 
    dqxi[7] = -dqxi[6]; 
    //---------------------
    std::vector<double> dqeta(8);
    dqeta[3] = (1.0 - xi) * (1.0 - zeta);
    dqeta[2] = xi * (1.0 - zeta);
    dqeta[6] = xi * zeta;
    dqeta[7] = (1.0 - xi) * zeta; 
    dqeta[0] = -dqeta[3];
    dqeta[1] = -dqeta[2];
    dqeta[5] = -dqeta[6];
    dqeta[4] = -dqeta[7];
    //---------------------
    std::vector<double> dqzeta(8);
    dqzeta[4] = (1.0 - xi) * (1.0 - eta);
    dqzeta[5] = xi * (1.0 - eta);
    dqzeta[6] = xi * eta;
    dqzeta[7] = (1.0 - xi) * eta;
    dqzeta[0] = -dqzeta[4];
    dqzeta[1] = -dqzeta[5];
    dqzeta[2] = -dqzeta[6];
    dqzeta[3] = -dqzeta[7];
    //---------------------
    for (int ii = 0; ii < 8; ++ii) {
      std::vector<double> c = d_vertexList[ii]->GetCoordinates();
      //---
      for (int jj = 0; jj < 3; ++jj) {
        B[jj] += c[jj] * dqxi[ii];
        B[jj + d_dim] += c[jj] * dqeta[ii];
        B[jj + 2*d_dim] += c[jj] * dqzeta[ii];
      }
      //---
    } // for (int ii = 0; ii < 8; ++ii)
  } // if (d_dim == 3)

  return true;

}


double Hexahedron::GetMeasure
() const
{

  //--- Get the quadrature rule
  std::vector<double> points;
  std::vector<double> weights;
  MathCXX::getQuadrature<double>("Gauss", 2, points, weights);

  //--- Map from [-1, 1] to [0, 1]
  for (int ii = 0; ii < 2; ++ii) {
    points[ii] = 0.5*(1 + points[ii]);
    weights[ii] *= 0.5;
  }

  std::vector<double> lCoord(3);
  std::vector<double> B(9);

  double meas = 0.0;
  for (int ixi = 0; ixi < 2; ++ixi) {
    lCoord[0] = points[ixi];
    for (int ieta = 0; ieta < 2; ++ieta) {
      lCoord[1] = points[ieta];
      for (int izeta = 0; izeta < 2; ++izeta) {
        lCoord[2] = points[izeta];
        GetMapJacobian(lCoord, B);
        double detB = 0.0;
        detB += B[0] * (B[4] * B[8] - B[5] * B[7]);
        detB -= B[1] * (B[3] * B[8] - B[5] * B[6]);
        detB += B[2] * (B[3] * B[7] - B[4] * B[6]);
        detB = (detB < 0.0) ? -detB : detB;
        meas += detB * weights[ixi] * weights[ieta] * weights[izeta];
      } // for (int izeta = 0; izeta < 2; ++izeta)
    } // for (int ieta = 0; ieta < 2; ++ieta)
  } // for (int ixi = 0; ixi < 2; ++ixi)

  return meas;

}


std::vector<GeometricElement*> Hexahedron::GetSideList
() const 
{

  // Create list of new pointers to the sides.
  std::vector<GeometricElement*> list(6);
  std::vector<Vertex*> tmpList(4);

  //--- Normal (N0, N1) x (N0, N2) is pointing inside the element.
  tmpList[0] = d_vertexList[4];
  tmpList[1] = d_vertexList[5];
  tmpList[2] = d_vertexList[6];
  tmpList[3] = d_vertexList[7];
  list[0] = new Quadrilateral(d_dim, tmpList, d_onBoundary);
  
  //--- Normal (N0, N1) x (N0, N2) is pointing inside the element.
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[4];
  tmpList[2] = d_vertexList[7];
  tmpList[3] = d_vertexList[3];
  list[1] = new Quadrilateral(d_dim, tmpList, d_onBoundary);

  //--- Normal (N0, N1) x (N0, N2) is pointing inside the element.
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  tmpList[2] = d_vertexList[5];
  tmpList[3] = d_vertexList[4];
  list[2] = new Quadrilateral(d_dim, tmpList, d_onBoundary);

  //--- Normal (N0, N1) x (N0, N2) is pointing outside the element.
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  tmpList[2] = d_vertexList[2];
  tmpList[3] = d_vertexList[3];
  list[3] = new Quadrilateral(d_dim, tmpList, d_onBoundary);

  //--- Normal (N0, N1) x (N0, N2) is pointing outside the element.
  tmpList[0] = d_vertexList[1];
  tmpList[1] = d_vertexList[5];
  tmpList[2] = d_vertexList[6];
  tmpList[3] = d_vertexList[2];
  list[4] = new Quadrilateral(d_dim, tmpList, d_onBoundary);

  //--- Normal (N0, N1) x (N0, N2) is pointing outside the element.
  tmpList[0] = d_vertexList[3];
  tmpList[1] = d_vertexList[2];
  tmpList[2] = d_vertexList[6];
  tmpList[3] = d_vertexList[7];
  list[5] = new Quadrilateral(d_dim, tmpList, d_onBoundary);

  return list;

}


bool Hexahedron::CheckMatch
(
  GeometricElement *ref
) const 
{

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Hexahedron *refHex = dynamic_cast<Hexahedron*>(ref);
  if (refHex == 0)
    return false;

  bool match = true;
  for (int i = 0; i < 8; ++i) {
    Vertex *iVertex = d_vertexList[i];
    bool belongsToRef = false;
    for (int j = 0; j < 8; ++j) {
      Vertex *jVertex = refHex->d_vertexList[j];
      if ((jVertex == iVertex) || ((*jVertex) == (*iVertex))) {
        belongsToRef = true;
        break;
      }
    } // for (int j = 0; j < 8; ++j)
    if (belongsToRef == false) {
      match = false;
      break;
    }
  } // for (int i = 0; i < 8; ++i)

  return match;

}


