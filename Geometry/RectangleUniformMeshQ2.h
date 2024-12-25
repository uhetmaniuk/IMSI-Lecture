#ifndef FECODE_GEOMETRY_RECTANGLE_UNIFORM_MESH_Q2_H
#define FECODE_GEOMETRY_RECTANGLE_UNIFORM_MESH_Q2_H

#include "Geometry/SimpleGeometry.h"


namespace FECode {

  //! A class generating a Q2 uniform mesh on a rectangle (9 nodes).

  /// ...
  /// ...
  class RectangleUniformMeshQ2: public SimpleGeometry {

    protected:

      unsigned int nX;
      double Lx;
      unsigned int nY;
      double Ly;

    public:

      /// Constructor
      RectangleUniformMeshQ2(unsigned int n_X, double L_x, unsigned int n_Y, double L_y);

      /// \brief Returns a partition of elements across processors.
      ///
      /// The routine creates a partition of elements.
      /// When numPart is 4^n, the number of elements in X-direction is a*2^n,
      /// and the number of elements in Y-direction is b*2^n, the partition
      /// is made manually.
      /// Else the routine calls METIS to partition the domain.
      /// When METIS is not linked to the code, a linear distribution of elements
      /// is done.
      ///
      /// \param[in] numPart Number of partitions to create.
      /// \param[in,out] cpuToCell Cell connectivity.
      int PartitionCells(int numPart, Connectivity &cpuToCell);

      /// Returns a list of global nodes on the boundary.
      std::vector<int> GetBoundaryNodes(location face) const; 

      /// Puts a Dirichlet condition on the boundary.
      int PutDirichletBC(std::vector<int> &hasDir) const;

      /// Returns a label describing the domain.
      const std::string Label() const;

      /// Returns a characteristic length along the X-axis for the domain (INLINE).
      double Length_X() const { return Lx; };

      /// Returns a characteristic length along the Y-axis for the domain (INLINE).
      double Length_Y() const { return Ly; };

      /// Returns the largest mesh size along the X-axis for the domain (INLINE).
      double MeshSize_X() const { return (nX <= 0) ? 0.0 : Lx/nX; };

      /// Returns the largest mesh size along the Y-axis for the domain (INLINE).
      double MeshSize_Y() const { return (nY <= 0) ? 0.0 : Ly/nY; };

      /// Outputs the geometry in vtk format.
      int Output_VTK(const std::string &fileName) const;

  };

}

#endif
