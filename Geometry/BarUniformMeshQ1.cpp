#include "Geometry/BarUniformMeshQ1.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Geometry/Connectivity.h"
#include "Utilities/Communicator.h"

using namespace FECode;
using namespace std;


BarUniformMeshQ1::BarUniformMeshQ1(unsigned int nEle, double L_X)
                 : SimpleGeometry(1, nEle, nEle+1),
                   Lx(L_X), nX(nEle)
                 {

  if (Lx <= 0.0) {
    cerr << "\n BarUniformMeshQ1: Wrong dimensions !!!\n\n";
    bool goodDomainLengths = false;
    assert(goodDomainLengths == true);
    return;
  }

  // Define global cells
  vector<int> list(2);
  for (int iEle = 0; iEle < nX; ++iEle) {
    list[0] = iEle;
    list[1] = iEle+1;
    Omega.SetCell(iEle, list);
  }

  // Define global points 
  if (nX > 0) {
    double hX = Lx/nX;
    vector<double> coord(1);
    for (unsigned int iEle = 0; iEle <= nX; ++iEle) {
      coord[0] = iEle*hX;
      Omega.SetPoint(iEle, coord);
    }
  }

}


std::vector<int> BarUniformMeshQ1::GetBoundaryNodes(location face) const {

  std::vector<int> list;

  switch (face) {
    case East:
      list.push_back(nX);
      break;
    case West:
      list.push_back(0);
      break;
    default:
      break;
  }

  return list;

}


int BarUniformMeshQ1::PutDirichletBC(vector<int> &hasDir) const {

  if (Omega.NumPoints() > 0) {
    hasDir.push_back(0);
    hasDir.push_back(Omega.NumPoints() - 1);
  }

  return 0;

}


const std::string BarUniformMeshQ1::Label() const {

  stringstream name;

  name << " >> Domain = [0, " << Lx << "]\n";
  name << " >> Orthogonal mesh uniform per direction with Q1 elements\n";
  name << endl;
  name << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
  name << endl;

  return name.str();

}


int BarUniformMeshQ1::Output_VTK(const string &fileName) const {

  stringstream name;
  name << fileName << ".geo";

  ofstream fout((name.str()).c_str());
  if (!fout) {
	  std::cerr << "\n\n !! Output geometric file " << name.str() << " opening failed !! \n\n\n";
    return -1;
  }

  fout << "# vtk DataFile Version 3.0" << endl;
  fout << fileName << endl;
  fout << "ASCII" << endl;
  fout << "DATASET STRUCTURED_POINTS" << endl;
  fout << "DIMENSIONS " << nX + 1 << " 1 1 " << endl;
  fout << "ORIGIN 0.0 0.0 0.0" << endl;
  fout << "SPACING " << Lx/nX << " 1.0 1.0 " << endl;

  fout.close();

  return 0;

}

