#include "Geometry/Quadrilateral.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "Geometry/Edge.h"
#include "Geometry/Vertex.h"
#include "Math/MathCXXCore.hpp"

using namespace FECode;


Quadrilateral::Quadrilateral
(
  short dim, 
  std::vector<Vertex*> &vertexList,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, QUADRILATERAL)
{

  // Test whether the vertex list is of size 4.
  if (vertexList.size() != 4) {
    std::cerr << "\n !!! A quadrilateral requires four distinct vertices !!!\n\n";
    assert(vertexList.size() == 4);
  }

  // Copy the vertex list.
  d_vertexList.assign(vertexList.begin(), vertexList.end());

}


Quadrilateral::~Quadrilateral() {
}


std::vector<double> Quadrilateral::GetCoordinatesList() const
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


bool Quadrilateral::GetMapJacobian
(
  const std::vector<double> &lCoord,
  std::vector<double> &B
) const 
{

  if (B.size() != 2*d_dim)
  {
    B.resize(2*d_dim);
  }

  if (d_dim == 3) {
    std::cerr << "\n !! Quadrilateral::GetMapJacobian >> ";
    std::cerr << "The Jacobian in 3D is not implemented !!\n\n";
    assert(0 > 1);
    return false;
  }  

  if (d_dim == 2) {
    //
    // Compute the Jacobian for a quadrilateral in two dimensions
    //
    std::vector<double> c0 = d_vertexList[0]->GetCoordinates();
    std::vector<double> c1 = d_vertexList[1]->GetCoordinates();
    std::vector<double> c2 = d_vertexList[2]->GetCoordinates();
    std::vector<double> c3 = d_vertexList[3]->GetCoordinates();

    double xi = lCoord[0];
    double eta = lCoord[1];

    //---------------------
    std::vector<double> dqxi(4);
    dqxi[0] = eta - 1.0;
    dqxi[1] = 1.0 - eta;
    dqxi[2] = eta;
    dqxi[3] = -eta;
    //---------------------
    std::vector<double> dqeta(4);
    dqeta[0] = xi - 1.0;
    dqeta[1] = -xi;
    dqeta[2] = xi;
    dqeta[3] = 1.0 - xi;
    //---------------------
    //--- Get the row [dx/dxi, dx/deta]
    B[0] = c0[0]*dqxi[0] + c1[0]*dqxi[1] + c2[0]*dqxi[2] + c3[0]*dqxi[3];
    B[2] = c0[0]*dqeta[0] + c1[0]*dqeta[1] + c2[0]*dqeta[2] + c3[0]*dqeta[3];
    //--- Get the row [dy/dxi, dy/deta]
    B[1] = c0[1]*dqxi[0] + c1[1]*dqxi[1] + c2[1]*dqxi[2] + c3[1]*dqxi[3];
    B[3] = c0[1]*dqeta[0] + c1[1]*dqeta[1] + c2[1]*dqeta[2] + c3[1]*dqeta[3];
  } // if (d_dim == 2)

  return true;

}


double Quadrilateral::GetMeasure
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

  std::vector<double> lCoord(2);
  std::vector<double> B(4);
  
  double meas = 0.0;
  for (int ixi = 0; ixi < 2; ++ixi) {

    lCoord[0] = points[ixi];
    for (int ieta = 0; ieta < 2; ++ieta) {
      lCoord[1] = points[ieta];
      GetMapJacobian(lCoord, B);
      double detB = B[0] * B[3] - B[1] * B[2];
      detB = (detB < 0.0) ? -detB : detB;
      meas += detB * weights[ixi] * weights[ieta];
    } // for (int ieta = 0; ieta < 2; ++ieta)

  } // for (int ixi = 0; ixi < 2; ++ixi)

  return meas;

}


std::vector<GeometricElement*> Quadrilateral::GetSideList
() const 
{

  // Create list of new pointers to the sides.
  std::vector<GeometricElement*> list(4);
  std::vector<Vertex*> tmpList(2);

  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  list[0] = new Edge(d_dim, tmpList, d_onBoundary);

  tmpList[0] = d_vertexList[1];
  tmpList[1] = d_vertexList[2];
  list[1] = new Edge(d_dim, tmpList, d_onBoundary);

  tmpList[0] = d_vertexList[2];
  tmpList[1] = d_vertexList[3];
  list[2] = new Edge(d_dim, tmpList, d_onBoundary);

  tmpList[0] = d_vertexList[3];
  tmpList[1] = d_vertexList[0];
  list[3] = new Edge(d_dim, tmpList, d_onBoundary);

  return list;

}


bool Quadrilateral::CheckMatch
(
  GeometricElement *ref
) const 
{

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Quadrilateral *refQuad = dynamic_cast<Quadrilateral*>(ref);
  if (refQuad == 0)
    return false;

  bool match = true;
  for (int i = 0; i < 4; ++i) {
    Vertex *iVertex = d_vertexList[i];
    bool belongsToRef = false;
    for (int j = 0; j < 4; ++j) {
      Vertex *jVertex = refQuad->d_vertexList[j];
      if ((jVertex == iVertex) || ((*jVertex) == (*iVertex))) {
        belongsToRef = true;
        break;
      }
    } // for (int j = 0; j < 4; ++j)
    if (belongsToRef == false) {
      match = false;
      break;
    }
  } // for (int i = 0; i < 4; ++i)

  return match;

}


