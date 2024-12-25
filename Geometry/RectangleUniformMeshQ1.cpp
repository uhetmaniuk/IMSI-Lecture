#include "Geometry/RectangleUniformMeshQ1.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Geometry/Connectivity.h"
#include "Geometry/Partition.h"

using namespace std;
using namespace FECode;


RectangleUniformMeshQ1::RectangleUniformMeshQ1(unsigned int n_X, double L_x,
                                               unsigned int n_Y, double L_y)
                       : SimpleGeometry(2, n_X*n_Y, (n_X+1)*(n_Y+1)),
                         nX(n_X), Lx(L_x), nY(n_Y), Ly(L_y)
                       {

  if ((Lx <= 0.0) || (Ly <= 0.0)) {
    cerr << "\n RectangleUniformMeshQ1: Wrong dimensions !!!\n\n";
    bool goodDomainLengths = false;
    assert(goodDomainLengths == true);
    return;
  }

  // Define global cells
  vector<int> list(4);
  for (int iY = 0; iY < nY; ++iY) {
    for (int iX = 0; iX < nX; ++iX) {
      int bottomLeft = iX + (nX+1)*iY;
      list[0] = bottomLeft;
      list[1] = bottomLeft + 1;
      list[2] = bottomLeft + nX + 2;
      list[3] = bottomLeft + nX + 1;
      Omega.SetCell(iX + iY*nX, list);
    }
  }

  // Define global points 
  if (nX*nY > 0) {
    double hX = Lx/nX;
    double hY = Ly/nY;
    vector<double> coord(2);
    for (unsigned int iY = 0; iY <= nY; ++iY) {
      for (unsigned int iX = 0; iX <= nX; ++iX) {
        coord[0] = iX*hX;
        coord[1] = iY*hY;
        Omega.SetPoint(iX + iY*(nX+1), coord);
      }
    }
  }

}


std::vector<int> RectangleUniformMeshQ1::GetBoundaryNodes(location face) const {

  std::vector<int> list;

  switch (face) {
    case East:
      for (int iY = 0; iY <= nY; ++iY)
        list.push_back(nX + iY*(nX+1));
      break;
    case North:
      for (int iX = 0; iX <= nX; ++iX)
        list.push_back(iX + (nX+1)*nY);
      break;
    case South:
      for (int iX = 0; iX <= nX; ++iX)
        list.push_back(iX);
      break;
    case West:
      for (int iY = 0; iY <= nY; ++iY)
        list.push_back(iY*(nX+1));
      break;
    default:
      break;
  }

  return list;

}


int RectangleUniformMeshQ1::PutDirichletBC(vector<int> &hasDir) const {

  for (unsigned int iY = 0; iY <= nY; ++iY)
    for (unsigned int iX = 0; iX <= nX; ++iX)
      if ((iX % nX)*(iY % nY) == 0)
        hasDir.push_back(iX + iY*(nX+1));

  return 0;

}


const std::string RectangleUniformMeshQ1::Label() const {

  stringstream name;

  name << " >> Domain = [0, " << Lx << "] x [0, " << Ly << "]\n";
  name << " >> Orthogonal mesh uniform per direction with Q1 elements\n";
  name << endl;
  name << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
  name << " Number of elements in [0, " << Ly << "] (Y-direction): " << nY << endl;
  name << endl;

  return name.str();

}


int RectangleUniformMeshQ1::PartitionCells(int numPart, Connectivity &cpuToCell) {

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
    cellToCell.Compose(cellToNode, nodeToCell, 2);

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


int RectangleUniformMeshQ1::Output_VTK(const string &fileName) const {

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
  fout << "DIMENSIONS " << nX + 1 << " " << nY + 1 << " 1 " << endl;
  fout << "ORIGIN 0.0 0.0 0.0" << endl;
  fout << "SPACING " << Lx/nX << " " << Ly/nY << " 1.0 " << endl;

  fout.close();

  return 0;

}

