#include "Geometry/Edge.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Geometry/Vertex.h"

using namespace FECode;


Edge::Edge
(
  short dim, 
  std::vector<Vertex*> &vertexList,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, EDGE)
{

  // Test whether the vertex list is of size 2.
  if (vertexList.size() != 2) {
    std::cerr << "\n !!! An edge requires two distinct vertices !!!\n\n";
    assert(vertexList.size() == 2);
  }

  // Copy the vertex list.
  d_vertexList.assign(vertexList.begin(), vertexList.end());

}


Edge::~Edge() {
}


std::vector<double> Edge::GetCoordinatesList() const
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


double Edge::GetMeasure() const 
{

  double length = 0.0;

  std::vector<double> n0 = d_vertexList[0]->GetCoordinates();
  std::vector<double> n1 = d_vertexList[1]->GetCoordinates();

  for (int iD = 0; iD < n0.size(); ++iD) {
    length += (n0[iD] - n1[iD]) * (n0[iD] - n1[iD]);
  }
  length = sqrt(length);

  return length;

}


std::vector<GeometricElement*> Edge::GetSideList
() const 
{

  // An edge is defined by a vertex list.
  // Here, we copy directly d_vertexList.

  std::vector<GeometricElement*> list(2);
  list[0] = d_vertexList[0];
  list[1] = d_vertexList[1];

  return list;

}


bool Edge::CheckMatch
(
  GeometricElement *ref
) const 
{

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Edge *refEdge = dynamic_cast<Edge*>(ref);
  if (refEdge == 0)
    return false;

  bool match = true;
  for (int i = 0; i < 2; ++i) {
    Vertex *iVertex = d_vertexList[i];
    bool belongsToRef = false;
    for (int j = 0; j < 2; ++j) {
      Vertex *jVertex = refEdge->d_vertexList[j];
      if ((jVertex == iVertex) || ((*jVertex) == (*iVertex))) {
        belongsToRef = true;
        break;
      }
    } // for (int j = 0; j < 2; ++j)
    if (belongsToRef == false) {
      match = false;
      break;
    }
  } // for (int i = 0; i < 2; ++i)

  return match;

}


