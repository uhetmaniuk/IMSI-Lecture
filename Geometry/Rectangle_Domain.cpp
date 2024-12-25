#include "Geometry/Rectangle_Domain.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


#include "Geometry/Mesh.h"
#include "Geometry/Vertex.h"
#include "Geometry/Quadrilateral.h"
#include "Geometry/Partition.h"

#include "Utilities/Communicator.h"
#include "Utilities/ParameterList.h"


using namespace FECode;


Rectangle_Domain::Rectangle_Domain
(
  const Communicator &Comm,
  const ParameterList &ParamPb
)
: 
  GeometricDomain(Comm, 2)
  , d_numEleX(0), d_numEleY(0)
  , d_ax(0.0), d_bx(0.0)
  , d_ay(0.0), d_by(0.0)
{

  //--- Read the parameters
  if (ParamPb.Get("Domain::Dimension_X", d_bx) == false) 
  {
    std::cerr << "\n !!! The dimension in x-direction is not specified !!! \n\n";
    assert(0 > 1);
  }

  if (ParamPb.Get("Domain::Dimension_Y", d_by) == false) 
  {
    std::cerr << "\n !!! The dimension in y-direction is not specified !!! \n\n";
    assert(0 > 1);
  }

  if ((d_bx > d_ax) && (d_by > d_ay))
  {
    std::stringstream name;
    name << " >> Domain = [" << d_ax << ", " << d_bx << "] ";
    name << " x [" << d_ay << ", " << d_by << "] ";
    name << "\n";
    d_domainLabel = name.str();  
  }
  else
  {
    if (d_bx <= d_ax)
      std::cerr << "\n !!! The dimension in x-direction must be positive !!! \n\n";
    if (d_by <= d_ay)
      std::cerr << "\n !!! The dimension in y-direction must be positive !!! \n\n";
    assert(d_bx > d_ax);
    assert(d_by > d_ay);
  }

  //--- Create the global mesh
  std::string meshType;
  ParamPb.Get("Domain::MeshType", meshType);
  if (meshType == "RECTANGULAR")
  {
    this->MakeRectangularMesh(ParamPb);
  }
  else
  {
    std::cerr << " Domain::MeshType " << meshType << "\n";
    std::cerr << "\n !!! There is no recognized mesh type for this domain !!!\n\n";
    assert(0 > 1);
  }

  //--- Partition the domain
  //--- Distribute the mesh among processors
  PartitionDomain();

}


Rectangle_Domain::~Rectangle_Domain
(
)
{

}


void Rectangle_Domain::MakeRectangularMesh
(
  const ParameterList &ParamPb
)
{

  //--- Read the parameters
  if (ParamPb.Get("Domain::N_Elements_X", d_numEleX) == false) 
  {
    std::cerr << "\n !!! The number of elements in x-direction is not specified !!! \n\n";
    assert(0 > 1);
  }

  if (ParamPb.Get("Domain::N_Elements_Y", d_numEleY) == false) 
  {
    std::cerr << "\n !!! The number of elements in y-direction is not specified !!! \n\n";
    assert(0 > 1);
  }

  int numGlobalElements = d_numEleX * d_numEleY;
  int numGlobalNodes = (d_numEleX + 1) * (d_numEleY + 1);

  if ((d_numEleX > 0) && (d_numEleY > 0))
  {
    std::stringstream name;
    name << d_domainLabel;
    name << " >> Orthogonal mesh uniform per direction with rectangles\n";
    name << " Number of elements in [0, " << d_bx << "] (X-direction): " << d_numEleX << std::endl;
    name << " Number of elements in [0, " << d_by << "] (Y-direction): " << d_numEleY << std::endl;
    name << " Number of global elements: " << numGlobalElements << std::endl;
    name << " Number of global elements: " << numGlobalNodes << std::endl;
    d_domainLabel = name.str();  
  }
  else
  {
    if (d_numEleX <= 0)
      std::cerr << "\n !!! The number of elements in x-direction must be positive !!! \n\n";
    assert(d_numEleX > 0);
    if (d_numEleY <= 0)
      std::cerr << "\n !!! The number of elements in y-direction must be positive !!! \n\n";
    assert(d_numEleY > 0);
  }

  // Define the global mesh
  d_Mesh_p = new Mesh(numGlobalNodes, numGlobalElements);

  // Define global nodes 
  double hX = (d_bx - d_ax) / d_numEleX;
  double hY = (d_by - d_ay) / d_numEleY;
  std::vector<double> coord(2);
  for (int iY = 0; iY <= d_numEleY; ++iY) {
    for (int iX = 0; iX <= d_numEleX; ++iX) {
      coord[0] = iX * hX;
      coord[1] = iY * hY;
      Vertex* thisVertex_p = new Vertex(2, coord, false);
      d_Mesh_p->InsertVertex(thisVertex_p);
    }
  }
  
  // Define the global elements
  std::vector<Vertex*> vertexList(4);
  for (int iY = 0; iY < d_numEleY; ++iY) {
    for (int iX = 0; iX < d_numEleX; ++iX) {
      int bottomLeft = iX + (d_numEleX + 1) * iY;
      vertexList[0] = d_Mesh_p->GetVertex(bottomLeft);
      vertexList[1] = d_Mesh_p->GetVertex(bottomLeft + 1);
      vertexList[2] = d_Mesh_p->GetVertex(bottomLeft + d_numEleX + 2);
      vertexList[3] = d_Mesh_p->GetVertex(bottomLeft + d_numEleX + 1);
      Quadrilateral *thisQuad_p = new Quadrilateral(2, vertexList);
      d_Mesh_p->InsertCell(thisQuad_p);
    } // for (int iX = 0; iX < numEleX; ++iX)
  } //  for (int iY = 0; iY < numEleY; ++iY)

}


void Rectangle_Domain::PartitionDomain
(
)
{

  int numGlobalElements = d_Mesh_p->NumberCells();
  std::vector<int> cellToPart(numGlobalElements, 0);
  int numPart = d_Comm.NumProc();

  //--- Check whether numPart is a power of 4 
  //--- and whether d_numEleX and d_numEleY are proportional to power of 2.
  double NumPartRoot4 = floor(0.5 + log((double) numPart)/log(4.0));
  double NumPartPower4 = floor(0.5 + pow(4.0, NumPartRoot4));
  int factor = (int) floor(0.5 + pow(2.0, NumPartRoot4));
  
  if ((numPart == (int) NumPartPower4) && 
      (d_numEleX % factor == 0) && (d_numEleY % factor == 0)) {
    // Partition the rectangle
    int ratioX = d_numEleX / factor;
    int ratioY = d_numEleY / factor;
    int count = 0;
    for (int iy = 0; iy < d_numEleY; ++iy)
      for (int ix = 0; ix < d_numEleX; ++ix)
        cellToPart[count++] = (ix / ratioX) + factor*(iy / ratioY);
  }
  else 
  {
    GeometricDomain::DefaultPartitionDomain(cellToPart);
  }

  //--- Distribute the cells and vertices among processors
  //--- Note that the global mesh will be destroyed.
  GeometricDomain::DistributeGlobalMesh(cellToPart);
  
}


