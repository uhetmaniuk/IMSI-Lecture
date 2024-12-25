#include "Geometry/SimpleGeometry.h"

#include "Geometry/Connectivity.h"
#include "Geometry/Partition.h"
#include "Elements/LagrangeVFE.h"

using namespace FECode;


SimpleGeometry::SimpleGeometry(unsigned char d, int numEle, int numNodes)
               : Omega(d, numEle, numNodes)
               {
}

int SimpleGeometry::PartitionCells(int numPart, Connectivity &cpuToCell) {

  int globalNumEle = Omega.NumCells();

  int *pointer = new int[globalNumEle+1];
  pointer[0] = 0;
  for (int iE = 0; iE < globalNumEle; ++iE)
    pointer[iE+1] = pointer[iE] + 1;

  int *cellToPart = new int[globalNumEle];

#ifdef FEC_USE_METIS
  const Connectivity &cellToNode = Omega.GetCellConnectivity();
  int degree = 0;
  if (degree == 0) {
    // Get the polynomial degree of the first local element.
    // It assumes that all the elements have the same degree.
    LagrangeVFE firstFE(Omega.GetCellPoints(0), Omega.SpaceDimension());
    degree = firstFE.PolynomialDegree();
  }

  Connectivity nodeToCell;
  nodeToCell.IsReverseOf(cellToNode);

  Connectivity cellToCell;
  cellToCell.Compose(cellToNode, nodeToCell, (Omega.SpaceDimension() - 1)*degree + 1);

  nodeToCell.Clear();   

  Partitioner tool;
  tool.METISCellPartition(numPart, globalNumEle, cellToPart, 
                          cellToCell.GetPointer(), cellToCell.GetAdjacency());
#else
  Partitioner tool;
  tool.UniformCellPartition(numPart, globalNumEle, cellToPart);
#endif

  // Define the connectivity cellToCpu
  Connectivity cellToCpu;
  cellToCpu.MakeView(globalNumEle, pointer, globalNumEle, cellToPart);

  cpuToCell.IsReverseOf(cellToCpu);

  delete[] pointer;
  delete[] cellToPart;

  return 0;

}

