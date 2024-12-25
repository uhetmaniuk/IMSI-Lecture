#ifndef FECODE_DOMAIN_H
#define FECODE_DOMAIN_H

#include <vector>
#include <ostream>


#include "Geometry/Connectivity.h"


namespace FECode {

  //! A class to manage the mesh on a domain.

  /// To be done
  /// To be done
  class Domain {

    private:

      unsigned char spaceDim;
      
      int numCell;
      int lMax;
      std::vector<int> *cellList;
      Connectivity *cellToNodes;
      
      std::vector<double> nodesList;

    public:

      /// Constructor.
      Domain(unsigned char d = 0, int numEle = 0, int numNodes = 0);

      /// Destructor.
      ~Domain();

      /// Sets the space dimension for the domain.
      void SetSpaceDimension(unsigned char d) { spaceDim = d; };

      /// Gets the space dimension of the domain.
      unsigned char SpaceDimension() const { return spaceDim; };

      /// Sets a cell in the mesh.
      void SetCell(const int iE, const std::vector<int> &newCell);

      /// Gets a cell from the mesh.
      std::vector<int> GetCell(const unsigned int i) const;

      /// Gets the nodal coordinates of points in a cell.
      std::vector<double> GetCellPoints(const int i) const;

      /// Gets the connectivity of a cell.
      const Connectivity & GetCellConnectivity();

      /// Returns the size of a cell.
      int GetCellSize(const unsigned int i) const;

      /// Sets the total number of cells.
      void SetNumCells(int cellSize);

      /// Gets the local number of cells.
      int NumCells() const { return numCell; };

      /// Adds a point to the mesh.
      int AddPoint(const std::vector<double> &newPoint);

      /// Gets the nodal coordinates for a point.
      std::vector<double> GetPoint(const unsigned int i) const
            { return std::vector<double>(&nodesList[spaceDim*i], 
                                         &nodesList[spaceDim*i] + spaceDim); };

      /// Gets the number of points in the mesh. 
      int NumPoints() const { return nodesList.size() / spaceDim; };

      /// Sets the number of points in the mesh.
      void SetNumPoints(int pointSize);

      /// Sets a point in the mesh.
      int SetPoint(const int iN, const std::vector<double> &newPoint);

      /// Clears the storage for the mesh.
      void Clear();

      /// Overloads the output operator <<.
      friend std::ostream & operator<<(std::ostream & os, const Domain & dom);

  };

  /// Overloads the output operator <<.
  std::ostream& operator<<(std::ostream & os, const FECode::Domain & dom);

}

#endif

