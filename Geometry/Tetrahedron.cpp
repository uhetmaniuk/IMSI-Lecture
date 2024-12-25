#include "Geometry/Tetrahedron.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "Geometry/Triangle.h"
#include "Geometry/Vertex.h"

using namespace FECode;


Tetrahedron::Tetrahedron
(
  short dim, 
  std::vector<Vertex*> &vertexList,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, TETRAHEDRON)
{

  // Test whether the dimension is at least 3
  if (dim < 3) {
    std::cerr << "\n !!! A tetrahedron requires three dimensions !!!\n\n";
    assert(dim >= 3);
  }

  // Test whether the vertex list is of size 4.
  if (vertexList.size() != 4) {
    std::cerr << "\n !!! A tetrahedron requires four distinct vertices !!!\n\n";
    assert(vertexList.size() == 4);
  }

  // Copy the vertex list.
  d_vertexList.assign(vertexList.begin(), vertexList.end());

}


Tetrahedron::~Tetrahedron() {
}


std::vector<double> Tetrahedron::GetCoordinatesList() const
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


bool Tetrahedron::GetMapJacobian
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
    std::cerr << " !! Tetrahedron::GetMapJacobian >> Not implemented for high dimensions !!";
    std::cerr << std::endl;
    return false;
  }

  // Note the Jacobian is constant for a tetrahedron.
  // So we do not use the local coordinates.

  std::vector<double> c0 = d_vertexList[0]->GetCoordinates();
  std::vector<double> c1 = d_vertexList[1]->GetCoordinates();
  std::vector<double> c2 = d_vertexList[2]->GetCoordinates();
  std::vector<double> c3 = d_vertexList[3]->GetCoordinates();

  for (int ii = 0; ii < d_dim; ++ii) {
    double c0_ii = c0[ii];
    B[ii] = c1[ii] - c0_ii;
    B[ii + d_dim] = c2[ii] - c0_ii;
    B[ii + 2*d_dim] = c3[ii] - c0_ii;
  } // for (int ii = 0; ii < d_dim; ++ii)

  return true;

}


double Tetrahedron::GetMeasure
() const
{

  std::vector<double> B(9);
  std::vector<double> coord(3);

  GetMapJacobian(coord, B);
  //double detB = B.Determinant();
  double detB = 0.0;
  detB += B[0] * (B[4] * B[8] - B[5] * B[7]);
  detB -= B[1] * (B[3] * B[8] - B[5] * B[6]);
  detB += B[2] * (B[3] * B[7] - B[4] * B[6]);
  detB = detB / 6.0;

  return ((detB < 0.0) ? -detB : detB );

}


std::vector<GeometricElement*> Tetrahedron::GetSideList
() const
{
  // Create list of new pointers to the sides.
  std::vector<GeometricElement*> list(4);
  std::vector<Vertex*> tmpList(3);

  //-- Normal (01) x (02) is pointing inside the element
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  tmpList[2] = d_vertexList[2];
  list[0] = new Triangle(d_dim, tmpList, d_onBoundary);

  //-- Normal (02) x (03) is pointing inside the element
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[2];
  tmpList[2] = d_vertexList[3];
  list[1] = new Triangle(d_dim, tmpList, d_onBoundary);

  //-- Normal (01) x (03) is pointing outside the element
  tmpList[0] = d_vertexList[0];
  tmpList[1] = d_vertexList[1];
  tmpList[2] = d_vertexList[3];
  list[2] = new Triangle(d_dim, tmpList, d_onBoundary);

  //-- Normal (12) x (13) is pointing outside the element
  tmpList[0] = d_vertexList[1];
  tmpList[1] = d_vertexList[2];
  tmpList[2] = d_vertexList[3];
  list[3] = new Triangle(d_dim, tmpList, d_onBoundary);

  return list;

}


bool Tetrahedron::CheckMatch
(
  GeometricElement *ref
) const 
{

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Tetrahedron *refTet = dynamic_cast<Tetrahedron*>(ref);
  if (refTet == 0)
    return false;

  bool match = true;
  for (int i = 0; i < 4; ++i) {
    Vertex *iVertex = d_vertexList[i];
    bool belongsToRef = false;
    for (int j = 0; j < 4; ++j) {
      Vertex *jVertex = refTet->d_vertexList[j];
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


