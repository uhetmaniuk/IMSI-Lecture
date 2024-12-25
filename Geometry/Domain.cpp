#include "Geometry/Domain.h"

#include <cstring>

#include "Geometry/Connectivity.h"

using namespace FECode;
using namespace std;


Domain::Domain(unsigned char d, int numEle, int numNodes)
       : spaceDim(d),
         numCell(numEle),
         lMax(0), cellList(0),
         cellToNodes(0), 
         nodesList(d*numNodes) {
         
  if (numCell > 0) {
    lMax = ((numCell + 7) / 8) * 8;
    cellList = new vector<int>[lMax];
  }

}


Domain::~Domain() {

  if (cellList)
    delete[] cellList;
  cellList = 0;
  lMax = 0;
  
  if (cellToNodes)
    delete cellToNodes;
  cellToNodes = 0;

}


void Domain::SetCell(const int iE, const std::vector<int> &newCell) {

  if (cellToNodes) {
    // Convert the storage back to cellList
    lMax = ((numCell + 7) / 8) * 8;
    cellList = new vector<int>[lMax];
    int *myPointer = cellToNodes->GetPointer();
    int *myTarget = cellToNodes->GetAdjacency();
    for (int iC = 0; iC < numCell; ++iC) {
      int cSize = myPointer[iC+1] - myPointer[iC];
      int *icTarget = myTarget + myPointer[iC];
      vector<int> tmp(cSize);
      for (int j = 0; j < cSize; ++j)
        tmp[j] = icTarget[j];
      cellList[iC] = tmp;
    }
    delete cellToNodes;
    cellToNodes = 0;
  }

  if (iE < lMax) {
    cellList[iE] = newCell;
    numCell = (iE < numCell) ? numCell : iE + 1;
  }
  else {
    // Need to resize the array
    numCell = iE + 1;
    int newLMAX = ((numCell + 7) / 8) * 8;
    vector<int> *tmp = new vector<int>[newLMAX];
    for (int iL = 0; iL < lMax; ++iL)
      tmp[iL] = cellList[iL];
    tmp[iE] = newCell;
    delete[] cellList;
    cellList = tmp;
    lMax = newLMAX;
  }

}


std::vector<int> Domain::GetCell(const unsigned int i) const {

  if (cellToNodes) {
    vector<int> tmp;
    int *myPointer = cellToNodes->GetPointer();
    int iSize = myPointer[i+1] - myPointer[i];
    tmp.resize(iSize);
    int *iTarget = cellToNodes->GetAdjacency() + myPointer[i];
    for (int j = 0; j < iSize; ++j)
      tmp[j] = iTarget[j];
    return tmp;
  }
  else {
    return cellList[i];
  }

}


const Connectivity& Domain::GetCellConnectivity() {

  if (cellToNodes)
    return *cellToNodes;

  // Count the total length of adjacency
  int totalAdj = 0;
  for (int iC = 0; iC < numCell; ++iC)
    totalAdj += cellList[iC].size();

  cellToNodes = new Connectivity(numCell, totalAdj);
  int *myPointer = cellToNodes->GetPointer();
  int *myTarget = cellToNodes->GetAdjacency();

  myPointer[0] = 0;
  for (int iC = 0; iC < numCell; ++iC) {
    const vector<int> &icList = cellList[iC];
    int cSize = icList.size();
    myPointer[iC+1] = myPointer[iC] + cSize;
    int *icTarget = myTarget + myPointer[iC];
    for (int j = 0; j < cSize; ++j)
      icTarget[j] = icList[j];
  }

  delete[] cellList;
  cellList = 0;
  lMax = 0;

  return *cellToNodes; 
 
}


vector<double> Domain::GetCellPoints(const int iC) const {

  vector<double> list;

  // Check if the data is stored in cellToNodes or in cellList
  if (cellToNodes) {
    vector<int> tmp;
    int *myPointer = cellToNodes->GetPointer();
    int iSize = myPointer[iC+1] - myPointer[iC];
    if (iSize == 0)
      return list;
    tmp.resize(iSize);
    int *iTarget = (*cellToNodes)[iC];
    for (int j = 0; j < iSize; ++j)
      tmp[j] = iTarget[j];
    list.resize(spaceDim*iSize);
    for (int j = 0; j < iSize; ++j) {
      int nj = tmp[j];
      for (int jj = 0; jj < spaceDim; ++jj)
        list[j*spaceDim + jj] = nodesList[nj*spaceDim + jj];    
    }
  }
  else {
    vector<int> points = cellList[iC];
    int sizeP = points.size();
    if (sizeP == 0)
      return list;
    list.resize(spaceDim*sizeP);
    for (int j = 0; j < sizeP; ++j) {
      int nj = points[j];
      for (int jj = 0; jj < spaceDim; ++jj)
        list[j*spaceDim + jj] = nodesList[nj*spaceDim + jj];    
    }
  }

  return list;

}


int Domain::GetCellSize(const unsigned int i) const {

  if (cellToNodes) { 
    int *myPointer = cellToNodes->GetPointer();
    int iSize = myPointer[i+1] - myPointer[i];
    return iSize;
  }
  else {
    return cellList[i].size();
  } 

}


void Domain::SetNumCells(int cellSize) {

  if (cellList)
    delete[] cellList;
  cellList = 0;
  lMax = 0;
  
  if (cellToNodes)
    delete cellToNodes;
  cellToNodes = 0;

  numCell = cellSize;
  lMax = ((numCell + 7) / 8) * 8;
  cellList = new vector<int>[lMax];

}


int Domain::AddPoint(const vector<double> &newPoint) {

  if (newPoint.size() != spaceDim)
    return -1;

  for (unsigned int j = 0; j < spaceDim; ++j)
    nodesList.push_back(newPoint[j]);

  return 0;

}


void Domain::SetNumPoints(int pointSize) {

  nodesList.clear();
  nodesList.resize(pointSize*spaceDim);

}


int Domain::SetPoint(const int iN, const vector<double> &newPoint) {

  if (iN < 0)
    return -1;
    
  if (newPoint.size() != spaceDim)
    return -2;

  if (iN*spaceDim >= nodesList.size()) {
    vector<double> tmp = nodesList;
    nodesList.resize((iN+1)*spaceDim);
    for (int i = 0; i < tmp.size(); ++i)
      nodesList[i] = tmp[i];
  }

  unsigned int offSet = iN*spaceDim;
  for (unsigned int j = 0; j < spaceDim; ++j)
    nodesList[offSet + j] = newPoint[j];
    
  return 0;

}


void Domain::Clear() {

  if (cellList)
    delete[] cellList;
  cellList = 0;
  lMax = 0;
  
  if (cellToNodes)
    delete cellToNodes;
  cellToNodes = 0;

  nodesList.clear();

}


std::ostream& FECode::operator<<(std::ostream & os, const FECode::Domain & dom) {

  os << " --- Domain of space dimension " << ((int) dom.spaceDim) << " --- " << endl;
  if (dom.cellToNodes) {
    os << *(dom.cellToNodes);
  }
  else {
    os << " --- Cell connectivity of size " << dom.numCell << " --- " << endl;
    for (int iC = 0; iC < dom.numCell; ++iC) {
      vector<int> list = dom.cellList[iC];
      os << " Cell " << iC << ": ";
      for (int jN = 0; jN < list.size(); ++jN)
        os << list[jN] << " ";
      os << endl;
    }
  }
  os << endl;
  os << " --- Coordinates of nodes --- " << endl;
  unsigned int sizeN = dom.nodesList.size()/dom.spaceDim;
  for (unsigned int i = 0; i < sizeN; ++i) {
    os << " Node " << i << ": ";
    for (unsigned int j = 0; j < dom.spaceDim; ++j)
      os << dom.nodesList[dom.spaceDim*i + j] << " ";
    os << endl;
  }
  os << endl;

  return os;

}


