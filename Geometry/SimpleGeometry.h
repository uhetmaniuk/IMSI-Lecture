#ifndef FECODE_GEOMETRY_SIMPLE_H
#define FECODE_GEOMETRY_SIMPLE_H

#include <string>
#include <vector>

#include "Geometry/Domain.h"


namespace FECode {

  //! An enumeration to locate faces of simple geometries.
  enum location{Back, East, Front, North, South, West};

  class Communicator;
  class Domain;

  //! An abstract class for simple geometries.

  class SimpleGeometry {

    protected:
      Domain Omega;

      /// Constructor
      SimpleGeometry(unsigned char d = 0, int numEle = 0, int numNodes = 0);

    public:

      /// Destructor.
      virtual ~SimpleGeometry() { };

      /// \brief Returns a partition of elements across processors.
      ///
      /// The routine creates a partition of elements by calling METIS.
      /// When METIS is not linked to the code, a linear distribution of elements
      /// is done.
      ///
      /// \param[in] numPart Number of partitions to create.
      /// \param[in,out] cpuToCell Cell connectivity.
      virtual int PartitionCells(int numPart, Connectivity &cpuToCell);

      /// Puts a Dirichlet condition on the boundary.
      virtual int PutDirichletBC(std::vector<int> &hasDir) const = 0;

      /// Clears the storage from the mesh.
      virtual void Clear() { Omega.Clear(); };

      /// Returns a list of global nodes on the boundary.
      virtual std::vector<int> GetBoundaryNodes(location face) const = 0;

      /// Returns a reference to describe the domain.
      virtual const Domain & GetDomain() const
            { return Omega; };

      /// Returns the connectivity of cells for the mesh.
      virtual const Connectivity & GetCellConnectivity()
            { return Omega.GetCellConnectivity(); };

      /// Returns the number of cells for the mesh.
      virtual const unsigned int NumCells() const
            { return Omega.NumCells(); }; 

      /// Returns the number of points for the mesh.
      virtual const unsigned int NumPoints() const 
            { return Omega.NumPoints(); }; 

      /// Returns a label describing the domain.
      virtual const std::string Label() const = 0;

      /// Returns a characteristic length along the X-axis for the domain.
      virtual double Length_X() const { return 0.0; };

      /// Returns a characteristic length along the Y-axis for the domain.
      virtual double Length_Y() const { return 0.0; };

      /// Returns a characteristic length along the Z-axis for the domain.
      virtual double Length_Z() const { return 0.0; };

      /// Returns the largest mesh size along the X-axis for the domain.
      virtual double MeshSize_X() const { return 0.0; };

      /// Returns the largest mesh size along the Y-axis for the domain.
      virtual double MeshSize_Y() const { return 0.0; };

      /// Returns the largest mesh size along the Z-axis for the domain.
      virtual double MeshSize_Z() const { return 0.0; };

      /// Outputs the geometry in vtk format.
      virtual int Output_VTK(const std::string &fileName) const { return -1; };
      
      /// Returns the space dimension of the domain.
      virtual unsigned char SpaceDimension() const
              { return Omega.SpaceDimension(); };

  };

}


#endif
