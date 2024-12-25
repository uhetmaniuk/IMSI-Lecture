#include "Geometry/RectangleUniformMeshQ2.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Geometry/Connectivity.h"
#include "Geometry/Partition.h"

using namespace std;
using namespace FECode;


RectangleUniformMeshQ2::RectangleUniformMeshQ2(unsigned int n_X, double L_x,
                                               unsigned int n_Y, double L_y)
                       : SimpleGeometry(2, n_X*n_Y, (2*n_X+1)*(2*n_Y+1)),
                         nX(n_X), Lx(L_x), nY(n_Y), Ly(L_y)
                       {

  if ((Lx <= 0.0) || (Ly <= 0.0)) {
    cerr << "\n RectangleUniformMeshQ2: Wrong dimensions !!!\n\n";
    bool goodDomainLengths = false;
    assert(goodDomainLengths == true);
    return;
  }

  // Define global cells
  vector<int> list(9);
  for (int iY = 0; iY < nY; ++iY) {
    for (int iX = 0; iX < nX; ++iX) {
      int bottomLeft = 2*iX + (2*nX+1)*2*iY;
      list[0] = bottomLeft;
      list[1] = bottomLeft + 2;
      list[2] = bottomLeft + 4*nX + 4;
      list[3] = bottomLeft + 4*nX + 2;
      list[4] = bottomLeft + 1;
      list[5] = bottomLeft + 2*nX + 3;
      list[6] = bottomLeft + 4*nX + 3;
      list[7] = bottomLeft + 2*nX + 1;
      list[8] = bottomLeft + 2*nX + 2;
      Omega.SetCell(iX + iY*nX, list);
    }
  }

  // Define global points 
  if (nX*nY > 0) {
    double hX = Lx/nX;
    double hY = Ly/nY;
    vector<double> coord(2);
    unsigned int twoNX = 2*nX;
    unsigned int twoNY = 2*nY;
    for (unsigned int iY = 0; iY <= twoNY; ++iY) {
      for (unsigned int iX = 0; iX <= twoNX; ++iX) {
        coord[0] = 0.5*iX*hX;
        coord[1] = 0.5*iY*hY;
        Omega.SetPoint(iX + iY*(twoNX+1), coord);
      }
    }
  }

}


std::vector<int> RectangleUniformMeshQ2::GetBoundaryNodes(location face) const {

  std::vector<int> list;

  switch (face) {
    case East:
      for (int iY = 0; iY <= 2*nY; ++iY)
        list.push_back(2*nX + iY*(2*nX+1));
      break;           
    case North:          
      for (int iX = 0; iX <= 2*nX; ++iX)
        list.push_back(iX + (2*nX+1)*2*nY);
      break; 
    case South: 
      for (int iX = 0; iX <= 2*nX; ++iX)
        list.push_back(iX);
      break;
    case West:
      for (int iY = 0; iY <= 2*nY; ++iY)
        list.push_back(iY*(2*nX+1));
      break;
    default:
      break;
  }

  return list;

}


int RectangleUniformMeshQ2::PutDirichletBC(vector<int> &hasDir) const {

  unsigned int twoNX = 2*nX;
  unsigned int twoNY = 2*nY;
  unsigned int count = 0;
  for (unsigned int iY = 0; iY <= twoNY; ++iY)
    for (unsigned int iX = 0; iX <= twoNX; ++iX)
      if ((iX % twoNX)*(iY % twoNY) == 0)
        hasDir.push_back(iX + iY*(twoNX+1));
        
  return 0;

}


const std::string RectangleUniformMeshQ2::Label() const {

  stringstream name;

  name << " >> Domain = [0, " << Lx << "] x [0, " << Ly << "]\n";
  name << " >> Orthogonal mesh uniform per direction with Q2 elements (9 nodes)\n";
  name << endl;
  name << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
  name << " Number of elements in [0, " << Ly << "] (Y-direction): " << nY << endl;
  name << endl;

  return name.str();

}


int RectangleUniformMeshQ2::PartitionCells(int numPart, Connectivity &cpuToCell) {

  int globalNumEle = Omega.NumCells();

  int *pointer = new int[2*globalNumEle+1];
  pointer[0] = 0;
  for (int iE = 0; iE < globalNumEle; ++iE)
    pointer[iE+1] = pointer[iE] + 1;

  int *cellToPart = pointer + globalNumEle + 1;

  //--- Check whether numPart is a power of 4 
  //--- and whether nX and nY are proportional to power of 2.
  double NumPartRoot4 = floor(0.5+log((double) numPart)/log(4.0));
  double NumPartPower4 = floor(0.5 + pow(4.0, NumPartRoot4));
  int factor = (int) floor(0.5 + pow(2.0, NumPartRoot4));

  if ((numPart == (int) NumPartPower4) && (nX % factor == 0) && (nY % factor == 0)) {
    // Partition the rectangle
    int ratioX = nX / factor;
    int ratioY = nY / factor;
    int count = 0;
    for (int iy = 0; iy < nY; ++iy)
      for (int ix = 0; ix < nX; ++ix)
        cellToPart[count++] = (ix / ratioX) + factor*(iy / ratioY);
  }
  else {
#ifdef FEC_USE_METIS
    const Connectivity &cellToNode = Omega.GetCellConnectivity();

    Connectivity nodeToCell;
    nodeToCell.IsReverseOf(cellToNode);

    Connectivity cellToCell;
    cellToCell.Compose(cellToNode, nodeToCell, 3);

    nodeToCell.Clear();   

    Partitioner tool;
    tool.METISCellPartition(numPart, globalNumEle, cellToPart, 
                            cellToCell.GetPointer(), cellToCell.GetAdjacency());
#else
    Partitioner tool;
    tool.UniformCellPartition(numPart, globalNumEle, cellToPart);
#endif
  }

  // Define the connectivity cellToCpu
  Connectivity cellToCpu;
  cellToCpu.MakeView(globalNumEle, pointer, globalNumEle, cellToPart);

  cpuToCell.IsReverseOf(cellToCpu);

  delete[] pointer;

  return 0;

}


int RectangleUniformMeshQ2::Output_VTK(const string &fileName) const {

  stringstream name;
  name << fileName << ".geo";

  ofstream fout((name.str()).c_str());
  if (!fout) {
    cerr << "\n\n !! Output geometric file " << name.str() << " opening failed !! \n\n\n";
    return -1;
  }

  fout << "# vtk DataFile Version 3.0" << endl;
  fout << fileName << endl;
  fout << "ASCII" << endl;
  fout << "DATASET STRUCTURED_POINTS" << endl;
  fout << "DIMENSIONS " << 2*nX + 1 << " " << 2*nY + 1 << " 1 " << endl;
  fout << "ORIGIN 0.0 0.0 0.0" << endl;
  fout << "SPACING " << 0.5*Lx/nX << " " << 0.5*Ly/nY << " 1.0 " << endl;

  fout.close();

  return 0;

}


