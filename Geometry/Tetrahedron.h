#ifndef FECODE_TETRAHEDRON_H
#define FECODE_TETRAHEDRON_H

#include <vector>

#include "Geometry/GeometricElement.h"

namespace FECode {

  class Vertex;

  //! A concrete class to define a tetrahedron in a mesh.

  /// The Tetrahedron class defines a tetrahedron with straight edges
  /// in a mesh on a domain.
  /// The class is derived from the virtual class GeometricElement.
  /// The reference tetrahedron is the unit right tetrahedron.
  /// The vertices are numbered as follows:
  /// Vertex 0 = (0, 0, 0), 
  /// Vertex 1 = (1, 0, 0),
  /// Vertex 2 = (0, 1, 0),
  /// Vertex 3 = (0, 0, 1).
  ///
  class Tetrahedron: public virtual GeometricElement {

    private:

      // Don't define these functions
      Tetrahedron(const Tetrahedron &ref);
      Tetrahedron & operator=(const Tetrahedron &ref);

    protected:

    public:

      //////////////////////
      // ACCESSOR
      //////////////////////


      /// Returns the list of coordinates for all the vertices.
      std::vector<double> GetCoordinatesList() const;


      /// Computes the Jacobian of the map from the reference element
      /// to the physical element at a local point. 
      /// The reference element is here the unit right tetrahedron.
      /// \param[in] lCoord: Vector of local coordinates of the point.
      /// \param[in,out] B: Jacobian matrix of the map at the local point.
      /// \return Boolean flag, set to TRUE when the computation is possible
      /// else set to FALSE.
      /// \note In 3D, the Jacobian matrix is 
      /// [dx/di, dx/eta, dx/dzeta; 
      ///  dy/dxi, dy/deta, dy/dzeta;
      ///  dz/dxi, dz/deta, dz/dzeta].
 	    /// \note The tetrahedron is affine-equivalent. So lCoord is not used.
      bool GetMapJacobian
      (
        const std::vector<double> &lCoord,
        std::vector<double> &B
      ) const;

      /// Returns the measure of the element.
      double GetMeasure() const;
  
      /// Returns the list of side elements.
      /// \note For the first two sides, the triangle (N1, N2, N3) is such that
      /// (N2 - N1) x (N3 - N1) is pointing inside the element.
      /// \note For the last two sides, the triangle (N1, N2, N3) is such that
      /// (N2 - N1) x (N3 - N1) is pointing outside the element.
      std::vector<GeometricElement*> GetSideList() const;

      /// Returns the number of vertices in a tetrahedron.
      int NumberVertices() const { return 4; };

      /// Returns the number of sides in a tetrahedron.
      int NumberSides() const { return 4; };

      ////////////
      // CREATOR
      ////////////

      /// Constructor.
      /// \param[in] dim - Space dimension (short).
      /// \param[in] vertexList - List of vertices forming the tetrahedron.
      /// \param[in] boundFlag - Boolean to flag whether the tetrahedron 
      ///  is on the boundary (Default value = 'false').
      Tetrahedron
      (
        short dim, 
        std::vector<Vertex*> &vertexList,
        bool boundFlag = false
      );

      /// Destructor.
      virtual ~Tetrahedron();

      ///////////////
      // MANIPULATOR
      ///////////////

      /// Checks whether the two elements are matching.
      /// \param[in] ref: Pointer to the element that neeeds to be compared.
      /// \return Boolean value whether the two elements match.
      bool CheckMatch(GeometricElement *ref) const; 

  };

}

#endif

