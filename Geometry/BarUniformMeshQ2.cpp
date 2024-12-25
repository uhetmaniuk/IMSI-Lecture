#include "Geometry/BarUniformMeshQ2.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Geometry/Connectivity.h"
#include "Utilities/Communicator.h"

using namespace FECode;
using namespace std;


BarUniformMeshQ2::BarUniformMeshQ2(unsigned int nEle, double L)
                 : SimpleGeometry(1, nEle, 2*nEle+1),
                   Lx(L), nX(nEle)
                 {

  if (Lx <= 0.0) {
    cerr << "\n BarUniformMeshQ2: Wrong dimensions !!!\n\n";
    bool goodDomainLengths = false;
    assert(goodDomainLengths == true);
    return;
  }

  // Define global cells
  vector<int> list(3);
  for (int iEle = 0; iEle < nX; ++iEle) {
    list[0] = 2*iEle;
    list[1] = 2*(iEle+1);
    list[2] = 2*iEle + 1;
    Omega.SetCell(iEle, list);
  }

  // Define global points 
  if (nX > 0) {
    double hX = Lx/nX;
    vector<double> coord(1);
    for (unsigned int iEle = 0; iEle < nX; ++iEle) {
      coord[0] = iEle*hX;
      Omega.SetPoint(2*iEle, coord);
      coord[0] = (iEle + 0.5)*hX;
      Omega.SetPoint(2*iEle+1, coord);
    }
    coord[0] = Lx;
    Omega.SetPoint(2*nX, coord);
  }

}


std::vector<int> BarUniformMeshQ2::GetBoundaryNodes(location face) const {

  std::vector<int> list;

  switch (face) {
    case East:
      list.push_back(2*nX);
      break;
    case West:
      list.push_back(0);
      break;
    default:
      break;
  }

  return list;

}


int BarUniformMeshQ2::PutDirichletBC(vector<int> &hasDir) const {

  if (Omega.NumPoints() > 0) {
    hasDir.push_back(0);
    hasDir.push_back(Omega.NumPoints() - 1);
  }

  return 0;

}


const std::string BarUniformMeshQ2::Label() const {

  stringstream name;

  name << " >> Domain = [0, " << Lx << "]\n";
  name << " >> Orthogonal mesh uniform per direction with Q2 elements\n";
  name << endl;
  name << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
  name << endl;

  return name.str();

}


int BarUniformMeshQ2::Output_VTK(const string &fileName) const {

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
  fout << "DIMENSIONS " << 2*nX + 1 << " 1 1 " << endl;
  fout << "ORIGIN 0.0 0.0 0.0" << endl;
  fout << "SPACING " << 0.5*Lx/nX << " 1.0 1.0 " << endl;

  fout.close();

  return 0;

}


