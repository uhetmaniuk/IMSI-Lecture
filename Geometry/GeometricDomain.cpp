#include "Geometry/GeometricDomain.h"


#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>


#ifdef FEC_USE_TRILINOS
#include "Epetra_Map.h"
#include "Epetra_Util.h"
#endif


#include "Geometry/Mesh.h"
#include "Geometry/Partition.h"

#include "Math/Graph.h"

#include "Utilities/Communicator.h"


using namespace FECode;


GeometricDomain::GeometricDomain
(
 const Communicator &Comm
 , short dim
) :
d_Comm(Comm)
, d_dim(dim)
, d_Mesh_p(0)
, d_Global_CellToCell_p(0)
, d_Global_VertexToVertex_p(0)
#ifdef FEC_USE_TRILINOS
, d_CellMap_p(0)
, d_VertexMap_p(0)
, d_VertexMap_1to1_p(0)
#endif
, d_domainLabel()
{
}


GeometricDomain::~GeometricDomain
(
)
{
  
  if (d_Mesh_p)
    delete d_Mesh_p;
  d_Mesh_p = 0;
  
  if (d_Global_CellToCell_p)
    delete d_Global_CellToCell_p;
  d_Global_CellToCell_p = 0;
  
  if (d_Global_VertexToVertex_p)
    delete d_Global_VertexToVertex_p;
  d_Global_VertexToVertex_p = 0;

#ifdef FEC_USE_TRILINOS
  if (d_CellMap_p)
    delete d_CellMap_p;
  d_CellMap_p = 0;
  
  if (d_VertexMap_p)
    delete d_VertexMap_p;
  d_VertexMap_p = 0;
  
  if (d_VertexMap_1to1_p)
    delete d_VertexMap_1to1_p;
  d_VertexMap_1to1_p = 0;
#endif
  
}


void GeometricDomain::DefaultPartitionDomain
(
 std::vector<int> &cellToPart
 ) const
{
  
  int numGlobalCells = d_Mesh_p->NumberCells();
  if (cellToPart.size() != numGlobalCells)
    cellToPart.resize(numGlobalCells);
  int numPart = d_Comm.NumProc();
  
  {
    Partitioner tool;
#ifdef FEC_USE_METIS
    MathCXX::Graph cellToCell;
    d_Mesh_p->MakeCellToCell(cellToCell);
    tool.METISCellPartition(numPart, numGlobalCells, &(cellToPart[0]), 
                            cellToCell.GetPointer(), cellToCell.GetAdjacency());
#else
    tool.UniformCellPartition(numPart, numGlobalCells, &(cellToPart[0]));
#endif
  }
  
}


void GeometricDomain::DistributeGlobalMesh
(
 const std::vector<int> &cellToPart
 )
{
  
  int numGlobalCells = d_Mesh_p->NumberCells();
  int numPart = d_Comm.NumProc();
  
  //--- Get the list of elements in this processor
  int countCell = 0;
  for (int ii = 0; ii < numGlobalCells; ++ii) 
  {
    if (cellToPart[ii] == d_Comm.MyPID())
      countCell += 1;
  }
  
  std::vector<int> localCellID(countCell, 0);
  countCell = 0;
  for (int ii = 0; ii < numGlobalCells; ++ii) 
  {
    if (cellToPart[ii] == d_Comm.MyPID()) 
    {
      localCellID[countCell] = ii;
      countCell += 1;
    }
  }
  
#ifdef FEC_USE_TRILINOS
  d_CellMap_p = new Epetra_Map(numGlobalCells, countCell, &(localCellID[0]), 0, *(d_Comm.getEpetraObject()));
  assert(d_CellMap_p->NumGlobalElements() == numGlobalCells);
#else
  std::cerr << "\n !!! GeometricDomain::DistributeGlobalMesh is not implemented without Trilinos !!!\n\n";
  assert(0 > 1);
#endif
  
  //--- Define the map for the vertices
  std::vector<int> globalVertexID;
  d_Global_CellToCell_p = new MathCXX::Graph();
  d_Global_VertexToVertex_p = new MathCXX::Graph();
  {
    
    //--- Build the global cell-to-vertex connectivity
    MathCXX::Graph My_CellToVertex;
    d_Mesh_p->MakeCellToVertex(My_CellToVertex);
    
    //--- Build the global vertex-to-cell connectivity
    MathCXX::Graph My_VertexToCell;
    My_VertexToCell.IsReverseOf(My_CellToVertex);
    
    //--- Build the global cell-to-cell connectivity
    d_Global_CellToCell_p->Compose(My_CellToVertex, My_VertexToCell);

    //--- (07/10) We might remove the vertex-to-vertex connectivity.
    //--- The code to recompute the vertex-to-vertex connectivity is commented below.
    //--- It would return an Epetra_CrsGraph object.
    d_Global_VertexToVertex_p->Compose(My_VertexToCell, My_CellToVertex);
    
    //--- Get the list of vertices on this processor
    std::set<int> tmpVertexID;
    for (int ii = 0; ii < numGlobalCells; ++ii) 
    {
      if (cellToPart[ii] == d_Comm.MyPID()) 
      {
        int numVertices = My_CellToVertex.NumConnected(ii);
        int *vertexList = My_CellToVertex[ii];
        for (int jj = 0; jj < numVertices; ++jj)
        {
          tmpVertexID.insert(vertexList[jj]);
        }
      }
    }
    globalVertexID = std::vector<int>(tmpVertexID.begin(), tmpVertexID.end());
  
  }
  
#ifdef FEC_USE_TRILINOS
  d_VertexMap_p = new Epetra_Map(-1, globalVertexID.size(), &(globalVertexID[0]), 
                                 0, *(d_Comm.getEpetraObject()));
  //--- Make the one-to-one map for the vertices
  Epetra_Util MyTool;
  d_VertexMap_1to1_p = new Epetra_Map(MyTool.Create_OneToOne_Map(*d_VertexMap_p));
  assert(d_VertexMap_1to1_p->NumGlobalElements() == d_Mesh_p->NumberVertices());
#else
  std::cerr << "\n !!! GeometricDomain::DistributeGlobalMesh is not implemented without Trilinos !!!\n\n";
  assert(0 > 1);
#endif
  
  if (numPart == 1)
  {
    return;
  }
  
  //
  //--- Create a local mesh on the subdomain.
  //

  Mesh *MyLocalMesh_p = new Mesh(globalVertexID.size(), localCellID.size());
  
  //--- Extract the sub-mesh by transferring the vertices and cells.
  for (int iV = 0; iV < globalVertexID.size(); ++iV)
  {
    MyLocalMesh_p->InsertVertex(d_Mesh_p->GetVertex(globalVertexID[iV]));
    d_Mesh_p->SetVertex((Vertex*) 0, globalVertexID[iV]);
  }

  for (int iV = 0; iV < localCellID.size(); ++iV)
  {
    MyLocalMesh_p->InsertCell(d_Mesh_p->GetCell(localCellID[iV]));
    d_Mesh_p->SetCell((GeometricElement*) 0, localCellID[iV]);
  }

  //--- Destroy the global mesh
  delete d_Mesh_p;
  
  //--- Store the pointer to the local mesh
  d_Mesh_p = MyLocalMesh_p;
  
}


void GeometricDomain::PartitionDomain
(
)
{

  std::vector<int> cellToPart(d_Mesh_p->NumberCells(), 0);
  GeometricDomain::DefaultPartitionDomain(cellToPart);
  
  //--- Distribute the cells and vertices among processors
  //--- Note that the global mesh will be destroyed.
  GeometricDomain::DistributeGlobalMesh(cellToPart);
  
}


/*

//
// This code is commented so far. 
//

#ifdef FEC_USE_TRILINOS
Epetra_CrsGraph* GeometricDomain::Get_Global_VertexToVertex
(
) const
{
 
  if (d_VertexMap_1to1_p == 0)
  {
    std::cout << "\n !!! The 1-to-1 map of vertices is not defined !!!\n\n";
    return 0;
  }

  Epetra_FECrsGraph *vTov = new Epetra_FECrsGraph(Copy, *d_DofsMap_1to1_p, 0, true);
  
  {
    
    //--- Build the global cell-to-vertex connectivity
    MathCXX::Graph My_CellToVertex;
    d_Mesh_p->MakeCellToVertex(My_CellToVertex);
    
    int numCells = My_CellToVertex.Length();
    for (int ic = 0; ic < numCells; ++ic)
    {
      std::set<int> tmpVertexID;
      int *vertexList = My_CellToVertex[ic];
      int numVertices = My_CellToVertex.NumConnected(ic);
      //--- Convert the list of vertices to global indices
      for (int iv = 0; iv < numVertices; ++iv)
      {
        int gid = d_VertexMap_p->GID(vertexList[iv]);
        assert(gid > -1);
        vertexList[iv] = gid;
      }
      //--- Insert the connectivity between vertices
      for (int iv = 0; iv < numVertices; ++iv)
      {
        int info = vTov->InsertGlobalIndices(vertexList[iv], numVertices, vertexList);
      }
    }

  }
  
  vTov->GlobalAssemble();
  
  return vTov;
  
}
#else
#endif
 
*/


MathCXX::Graph* GeometricDomain::Get_Global_VertexToVertex
(
) const
{
  return d_Global_VertexToVertex_p;
}


const std::string GeometricDomain::Label
(
) const 
{
  
  std::stringstream name;
  
  name << d_domainLabel;
  name << std::endl;
  
  return name.str();
  
}


int GeometricDomain::NumberGlobalVertices
(
) const
{

#ifdef FEC_USE_TRILINOS
  
  return d_VertexMap_1to1_p->NumGlobalPoints();
  
#else

  return -1;

#endif

}


