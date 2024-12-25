#ifndef FECODE_GEOMETRY_BAR_UNIFORM_MESH_Q2_H
#define FECODE_GEOMETRY_BAR_UNIFORM_MESH_Q2_H

#include "Geometry/SimpleGeometry.h"

namespace FECode {

  //! A class generating a Q2 mesh on an interval.

  /// ...
  /// ...
  class BarUniformMeshQ2: public SimpleGeometry {

    protected:

      double Lx;
      unsigned int nX;

    public:

      /// Constructor
      BarUniformMeshQ2(unsigned int nEle, double L = 1.0);

      /// Returns a list of global nodes on the boundary.
      std::vector<int> GetBoundaryNodes(location face) const; 

      /// Puts a Dirichlet condition on the boundary.
      int PutDirichletBC(std::vector<int> &hasDir) const;

      /// Returns a label describing the domain.
      const std::string Label() const;

      /// Returns a characteristic length along the X-axis for the domain (INLINE).
      double Length_X() const { return Lx; };

      /// Returns the largest mesh size along the X-axis for the domain (INLINE).
      double MeshSize_X() const { return (nX <= 0) ? 0.0 : Lx/nX; };

      /// Outputs the geometry in vtk format.
      int Output_VTK(const std::string &fileName) const;

  };

}

#endif
