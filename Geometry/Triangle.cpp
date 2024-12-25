#include "Geometry/Triangle.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "Geometry/Edge.h"
#include "Geometry/Vertex.h"

using namespace FECode;


Triangle::Triangle
(
  short dim, 
  std::vector<Vertex*> &vertexList,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, TRIANGLE)
{

  // Test whether the vertex list is of size 3.
  if (vertexList.size() != 3) {
    std::cerr << "\n !!! A triangle requires three distinct vertices !!!\n\n";
    assert(vertexList.size() == 3);
  }

  // Copy the vertex list.
  d_vertexList.assign(vertexList.begin(), vertexList.end());

}


Triangle::~Triangle() {
}


std::vector<double> Triangle::GetCoordinatesList() const
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


bool Triangle::GetMapJacobian
(
  const std::vector<double> &lCoord, 
  std::vector<double> &B
) const 
{

  if (B.size() != 2*d_dim) 
  {
    B.resize(2*d_dim);
  }

  // Note the Jacobian is constant for a triangle.
  // So we do not use the local coordinates.

  std::vector<double> c0 = d_vertexList[0]->GetCoordinates();
  std::vector<double> c1 = d_vertexList[1]->GetCoordinates();
  std::vector<double> c2 = d_vertexList[2]->GetCoordinates();

  for (int ii = 0; ii < d_dim; ++ii) {
    double c0_ii = c0[ii];
    B[ii] = c1[ii] - c0_ii;
    B[ii + d_dim] = c2[ii] - c0_ii;
  } // for (int ii = 0; ii < d_dim; ++ii)

  return true;

}


double Triangle::GetMeasure() const
{

  std::vector<double> B(4);
  std::vector<double> coord(2);

  GetMapJacobian(coord, B);
  double detB = B[0] * B[3] - B[1] * B[2];

  return ((detB < 0.0) ? - 0.5 * detB : 0.5 * detB );

}


std::vector<GeometricElement*> Triangle::GetSideList
() const
{
  // Create list of new pointers to the sides.
  std::vector<GeometricElement*> list(3);
  std::vector<Vertex*> tmpList(2);

  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  list[0] = new Edge(d_dim, tmpList, d_onBoundary);

  tmpList[0] = d_vertexList[1];
  tmpList[1] = d_vertexList[2];
  list[1] = new Edge(d_dim, tmpList, d_onBoundary);

  tmpList[0] = d_vertexList[2];
  tmpList[1] = d_vertexList[0];
  list[2] = new Edge(d_dim, tmpList, d_onBoundary);

  return list;

}


bool Triangle::CheckMatch
(
  GeometricElement *ref
) const 
{

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Triangle *refTri = dynamic_cast<Triangle*>(ref);
  if (refTri == 0)
    return false;

  bool match = true;
  for (int i = 0; i < 3; ++i) {
    Vertex *iVertex = d_vertexList[i];
    bool belongsToRef = false;
    for (int j = 0; j < 3; ++j) {
      Vertex *jVertex = refTri->d_vertexList[j];
      if ((jVertex == iVertex) || ((*jVertex) == (*iVertex))) {
        belongsToRef = true;
        break;
      }
    } // for (int j = 0; j < 3; ++j)
    if (belongsToRef == false) {
      match = false;
      break;
    }
  } // for (int i = 0; i < 3; ++i)

  return match;

}


