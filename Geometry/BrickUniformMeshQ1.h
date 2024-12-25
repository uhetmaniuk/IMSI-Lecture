#ifndef FECODE_GEOMETRY_BRICK_UNIFORM_MESH_Q1_H
#define FECODE_GEOMETRY_BRICK_UNIFORM_MESH_Q1_H

#include "Geometry/SimpleGeometry.h"

namespace FECode {

  //! A class generating a Q1 uniform mesh on a brick.

  /// ...
  /// ...
  class BrickUniformMeshQ1: public SimpleGeometry {

    protected:

      unsigned int nX;
      double Lx;
      unsigned int nY;
      double Ly;
      unsigned int nZ;
      double Lz;

    public:

      /// Constructor
      BrickUniformMeshQ1(unsigned int n_X, double L_x, unsigned int n_Y, double L_y,
                         unsigned int n_Z, double L_z);

      /// \brief Returns a partition of elements across processors.
      ///
      /// The routine creates a partition of elements.
      /// When numPart is 8^n and the number of elements is a*2^n in X-direction,
      /// b*2^n in Y-direction, and c*2^n in Z-direction, the partition
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

      /// Returns a characteristic length along the Z-axis for the domain (INLINE).
      double Length_Z() const { return Lz; };

      /// Returns the largest mesh size along the X-axis for the domain (INLINE).
      double MeshSize_X() const { return (nX <= 0) ? 0.0 : Lx/nX; };

      /// Returns the largest mesh size along the Y-axis for the domain (INLINE).
      double MeshSize_Y() const { return (nY <= 0) ? 0.0 : Ly/nY; };

      /// Returns the largest mesh size along the Z-axis for the domain (INLINE).
      double MeshSize_Z() const { return (nZ <= 0) ? 0.0 : Lz/nZ; };

      /// Outputs the geometry in vtk format.
      int Output_VTK(const std::string &fileName) const;

  };

}

#endif
