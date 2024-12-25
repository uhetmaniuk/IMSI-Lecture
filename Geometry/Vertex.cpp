#include "Geometry/Vertex.h"

#include <vector>

using namespace FECode;


Vertex::Vertex
(
  short dim,
  std::vector<double> &coord,
  bool boundFlag
) :
  GeometricElement(dim, boundFlag, VERTEX),
  d_coordinates(dim, 0.0)
{

  // Copy the coordinates
  int nValues = (d_dim < coord.size()) ? (int) d_dim : coord.size();
  for (int ii = 0; ii < nValues; ++ii)
    d_coordinates[ii] = coord[ii];
    
}


Vertex::Vertex
(
  const Vertex &ref
) :
  GeometricElement(ref.d_dim, ref.d_onBoundary, VERTEX),
  d_coordinates(ref.d_dim, 0.0)
{

  // Copy the coordinates
  for (int ii = 0; ii < ref.d_dim; ++ii)
    d_coordinates[ii] = ref.d_coordinates[ii];
    
}


Vertex::~Vertex() {
}


std::vector<double> Vertex::GetCoordinates
() const 
{

  int size = d_coordinates.size();
  std::vector<double> myCoord(size, 0.0);
  for (int ii = 0; ii < size; ++ii)
    myCoord[ii] = d_coordinates[ii];

  return myCoord;

}


std::vector<Vertex*> Vertex::GetVertexList
() const 
{

  std::vector<Vertex*> list(1);
  list[0] = const_cast<Vertex*>(this);

  return list;

}


std::vector<GeometricElement*> Vertex::GetSideList
() const 
{

  std::vector<GeometricElement*> emptyList;
  return emptyList;

}


bool Vertex::CheckMatch(GeometricElement *ref) const {

  if (ref == 0)
    return false;

  if (this->d_Type != ref->GetType())
    return false;

  Vertex *refVertex = dynamic_cast<Vertex*>(ref);

  return ((*this) == (*refVertex));

}


int Vertex::SetCoordinates
(
  const std::vector<double> &coord
) 
{

  for (int iD = 0; iD < d_dim; ++iD)
    d_coordinates[iD] = coord[iD];

  return 0;

}


namespace FECode {
  
bool operator==
(
  const FECode::Vertex &n1, 
  const FECode::Vertex &n2
) 
{

  if (n1.d_dim != n2.d_dim)
    return false;

  bool match = true;
  for (int id = 0; id < n1.d_dim; ++id) {
    if (n1.d_coordinates[id] != n2.d_coordinates[id]) {
      match = false;
      break;
    }
  }

  return match;

}

}


