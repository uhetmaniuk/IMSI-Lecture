#include "Geometry/BrickUniformMeshQ2.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Geometry/Connectivity.h"
#include "Geometry/Partition.h"

using namespace std;
using namespace FECode;


BrickUniformMeshQ2::BrickUniformMeshQ2(unsigned int n_X, double L_x,
                                       unsigned int n_Y, double L_y,
                                       unsigned int n_Z, double L_z)
                   : SimpleGeometry(3, n_X*n_Y*n_Z, (2*n_X+1)*(2*n_Y+1)*(2*n_Z+1)),
                     nX(n_X), Lx(L_x), nY(n_Y), Ly(L_y), nZ(n_Z), Lz(L_z)
                   {

  if (nX*nY*nZ == 0) {
    cerr << "\n BrickUniformMeshQ2: Wrong number of elements !!!\n\n";
    bool goodNumberElements = false;
    assert(goodNumberElements == true);
    return;
  }
  
  if ((Lx <= 0.0) || (Ly <= 0.0) || (Lz <= 0.0)) {
    cerr << "\n BrickUniformMeshQ2: Wrong dimensions !!!\n\n";
    bool goodDomainLengths = false;
    assert(goodDomainLengths == true);
    return;
  }

  // Define global cells
  vector<int> list(27);
  for (int iZ = 0; iZ < nZ; ++iZ) {
    for (int iY = 0; iY < nY; ++iY) {
      for (int iX = 0; iX < nX; ++iX) {
        int bottomLeft = 2*iX + (2*nX+1)*2*iY + (2*nX+1)*(2*nY+1)*2*iZ;
        int middleLeft = bottomLeft + (2*nX+1)*(2*nY+1);
        int topLeft = middleLeft + (2*nX+1)*(2*nY+1);
        //----------------
        list[ 0] = bottomLeft;
        list[ 1] = bottomLeft + 2;
        list[ 2] = bottomLeft + 4*nX + 4;
        list[ 3] = bottomLeft + 4*nX + 2;
        //----------------
        list[ 4] = topLeft;
        list[ 5] = topLeft + 2;
        list[ 6] = topLeft + 4*nX + 4;
        list[ 7] = topLeft + 4*nX + 2;
        //----------------
        list[ 8] = bottomLeft + 1;
        list[ 9] = bottomLeft + 2*nX + 3;
        list[10] = bottomLeft + 4*nX + 3;
        list[11] = bottomLeft + 2*nX + 1;
        //----------------
        list[12] = topLeft + 1;
        list[13] = topLeft + 2*nX + 3;
        list[14] = topLeft + 4*nX + 3;
        list[15] = topLeft + 2*nX + 1;
        //----------------
        list[16] = middleLeft;
        list[17] = middleLeft + 2;
        list[18] = middleLeft + 4*nX + 4;
        list[19] = middleLeft + 4*nX + 2;
        //----------------
        list[20] = bottomLeft + 2*nX + 2;
        //----------------
        list[21] = topLeft + 2*nX + 2;
        //----------------
        list[22] = middleLeft + 1;
        //----------------
        list[23] = middleLeft + 4*nX + 3;
        //----------------
        list[24] = middleLeft + 2*nX + 1;
        //----------------
        list[25] = middleLeft + 2*nX + 3;
        //----------------
        list[26] = middleLeft + 2*nX + 2;
        //----------------
        Omega.SetCell(iX + iY*nX + iZ*nX*nY, list);
      }
    }
  }

  // Define global points
  if (nX*nY*nZ > 0) { 
    double hX = Lx/nX;
    double hY = Ly/nY;
    double hZ = Lz/nZ;
    int twoNX = 2*nX;
    int twoNY = 2*nY;
    int twoNZ = 2*nZ;
    vector<double> coord(3);
    for (int iZ = 0; iZ <= twoNZ; ++iZ) {
      for (int iY = 0; iY <= twoNY; ++iY) {
        for (int iX = 0; iX <= twoNX; ++iX) {
          coord[0] = 0.5*iX*hX;
          coord[1] = 0.5*iY*hY;
          coord[2] = 0.5*iZ*hZ;
          Omega.SetPoint(iX + iY*(twoNX+1) + iZ*(twoNX+1)*(twoNY+1), coord);
        }
      }
    }
  }

}


std::vector<int> BrickUniformMeshQ2::GetBoundaryNodes(location face) const {

  std::vector<int> list;

  int twoNX = 2*nX;
  int twoNY = 2*nY;
  int twoNZ = 2*nZ;

  switch (face) {
    case Back:
      for (int iY = 0; iY <= twoNY; ++iY)
        for (int iX = 0; iX <= twoNX; ++iX)
          list.push_back(iX + iY * (twoNX + 1));
      break;
    case East:
      for (int iZ = 0; iZ <= twoNZ; ++iZ)
        for (int iY = 0; iY <= twoNY; ++iY)
          list.push_back(twoNX + iY * (twoNX + 1) + iZ * (twoNX + 1) * (twoNY + 1));
      break;           
    case Front:
      for (int iY = 0; iY <= twoNY; ++iY)
        for (int iX = 0; iX <= twoNX; ++iX)
          list.push_back(iX + iY * (twoNX + 1) + (twoNX + 1) * (twoNY + 1) * twoNZ);
      break;
    case North:          
      for (int iZ = 0; iZ <= twoNZ; ++iZ)
        for (int iX = 0; iX <= twoNX; ++iX)
          list.push_back(iX + iZ * (twoNX + 1) * (twoNY + 1) + (twoNX + 1) * twoNY);
      break; 
    case South: 
      for (int iZ = 0; iZ <= twoNZ; ++iZ)
        for (int iX = 0; iX <= twoNX; ++iX)
          list.push_back(iX + iZ * (twoNX + 1) * (twoNY + 1));
      break;
    case West:
      for (int iZ = 0; iZ <= twoNZ; ++iZ)
        for (int iY = 0; iY <= twoNY; ++iY)
          list.push_back(iY * (twoNX + 1) + iZ * (twoNX + 1) * (twoNY + 1));
      break;
    default:
      break;
  }

  return list;

}


int BrickUniformMeshQ2::PutDirichletBC(vector<int> &hasDir) const {

  unsigned int twoNX = 2*nX;
  unsigned int twoNY = 2*nY;
  unsigned int twoNZ = 2*nZ;
  for (unsigned int iZ = 0; iZ <= twoNZ; ++iZ)
    for (unsigned int iY = 0; iY <= twoNY; ++iY)
      for (unsigned int iX = 0; iX <= twoNX; ++iX)
        if ((iX % twoNX)*(iY % twoNY)*(iZ % twoNZ) == 0)
          hasDir.push_back(iX + iY*(twoNX+1) + iZ*(twoNX+1)*(twoNY+1));

  return 0;

}


const std::string BrickUniformMeshQ2::Label() const {

  stringstream name;
  
  name << " >> Domain = [0, " << Lx << "] x [0, " << Ly << "] x [0, " << Lz << "]\n";
  name << " >> Orthogonal mesh uniform per direction with Q2 elements (27 nodes)\n";
  name << endl;
  name << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
  name << " Number of elements in [0, " << Ly << "] (Y-direction): " << nY << endl;
  name << " Number of elements in [0, " << Lz << "] (Z-direction): " << nZ << endl;
  name << endl;

  return name.str();

}


int BrickUniformMeshQ2::PartitionCells(int numPart, Connectivity &cpuToCell) {

  int globalNumEle = Omega.NumCells();

  int *pointer = new int[2*globalNumEle+1];
  pointer[0] = 0;
  for (int iE = 0; iE < globalNumEle; ++iE)
    pointer[iE+1] = pointer[iE] + 1;

  int *cellToPart = pointer + globalNumEle + 1;

  //--- Check whether numPart is a power of 8 
  //--- and whether nX, nY, and nZ are proportional to power of 2.
  double NumPartRoot8 = floor(0.5+log((double) numPart)/log(8.0));
  double NumPartPower8 = floor(0.5 + pow(8.0, NumPartRoot8));
  int factor = (int) floor(0.5 + pow(2.0, NumPartRoot8));

  if ((numPart == (int) NumPartPower8) && (nX % factor == 0) 
      && (nY % factor == 0) && (nZ % factor == 0)) {
    // Partition the rectangle
    int ratioX = nX / factor;
    int ratioY = nY / factor;
    int ratioZ = nZ / factor;
    int count = 0;
    for (int iz = 0; iz < nZ; ++iz)
      for (int iy = 0; iy < nY; ++iy)
        for (int ix = 0; ix < nX; ++ix)
          cellToPart[count++] = (ix / ratioX) + factor*(iy / ratioY) 
                                              + factor*factor*(iz / ratioZ);
  }
  else {
#ifdef FEC_USE_METIS
    const Connectivity &cellToNode = Omega.GetCellConnectivity();

    Connectivity nodeToCell;
    nodeToCell.IsReverseOf(cellToNode);

    Connectivity cellToCell;
    cellToCell.Compose(cellToNode, nodeToCell, 5);

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


int BrickUniformMeshQ2::Output_VTK(const string &fileName) const {

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
  fout << "DIMENSIONS " << 2*nX + 1 << " " << 2*nY + 1 << " " << 2*nZ + 1 << endl;
  fout << "ORIGIN 0.0 0.0 0.0" << endl;
  fout << "SPACING " << 0.5*Lx/nX << " " << 0.5*Ly/nY << " " << 0.5*Lz/nZ << endl;

  fout.close();

  return 0;

}


