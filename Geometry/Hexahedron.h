#ifndef FECODE_HEXAHEDRON_H
#define FECODE_HEXAHEDRON_H

#include <vector>

#include "Geometry/GeometricElement.h"

namespace FECode {
  
  class Vertex;
  
  //! A concrete class to define a hexahedron in a mesh.
  
  /// The Hexahedron class defines a hexahedron with straight edges
  /// in a mesh on a domain.
  /// The class is derived from the virtual class GeometricElement.
  /// The reference hexahedron is the unit cube [0, 1] x [0, 1] x [0, 1].
  /// The vertices are numbered as follows:
  /// Vertex 0 = (0, 0, 0),
  /// Vertex 1 = (1, 0, 0),
  /// Vertex 2 = (1, 1, 0),
  /// Vertex 3 = (0, 1, 0),
  /// Vertex 4 = (0, 0, 1),
  /// Vertex 5 = (1, 0, 1),
  /// Vertex 6 = (1, 1, 1),
  /// Vertex 7 = (0, 1, 1).
  ///
  class Hexahedron: public virtual GeometricElement {
    
  protected:
    
  public:
    
    //////////////////////
    // ACCESSOR
    //////////////////////


    /// Returns the list of coordinates for all the vertices.
    std::vector<double> GetCoordinatesList() const;


    /// Computes the Jacobian of the map from the reference element
    /// to the physical element at a local point.
    /// The reference element is here the unit cube [0,1] x [0,1] x [0,1].
    /// \param[in] lCoord: Vector of local coordinates of the point in [0,1] x [0,1] x [0,1].
    /// \param[out] B: Jacobian matrix of the map at the local point.
    /// \return Boolean flag, set to TRUE when the computation is possible
    /// else set to FALSE.
    /// \note In 3D, the Jacobian matrix is 
    /// [dx/di, dx/eta, dx/dzeta; 
    ///  dy/dxi, dy/deta, dy/dzeta;
    ///  dz/dxi, dz/deta, dz/dzeta].
    bool GetMapJacobian
    (
     const std::vector<double> &lCoord,
     std::vector<double> &B
     ) const;
    
    /// Returns the measure of the element.
    double GetMeasure() const;
    
    /// Returns the list of side elements.
    /// \note For the first three sides, the quadrilateral (N0, N1, N2, N3) is such that
    /// (N1 - N0) x (N2 - N0) is pointing inside the element.
    /// \note For the last three sides, the quadrilateral (N0, N1, N2, N3) is such that
    /// (N1 - N0) x (N2 - N0) is pointing outside the element.
    /// \note The normal outer to the hexahedron can vary
    /// on each quadrilateral side.
    std::vector<GeometricElement*> GetSideList() const;
    
    /// Returns the number of vertices in a hexahedron.
    int NumberVertices() const { return 8; };
    
    /// Returns the number of sides in a hexahedron. 
    int NumberSides() const { return 6; };
    
    ////////////
    // CREATOR
    ////////////
    
    /// Constructor.
    /// \param[in] dim - Space dimension (short).
    /// \param[in] vertexList - List of vertices forming the hexahedron.
    /// \param[in] boundFlag - Boolean to flag whether the hexahedron is on the boundary (Default value = 'false').
    Hexahedron
    (
     short dim, 
     std::vector<Vertex*> &vertexList,
     bool boundFlag = false
     );
    
    /// Destructor.
    virtual ~Hexahedron();
    
    ///////////////
    // MANIPULATOR
    ///////////////
    
    /// Checks whether the two elements are matching.
    /// \param[in] ref: Pointer to the element that neeeds to be compared.
    /// \return Boolean value whether the two elements match.
    bool CheckMatch
    (
     GeometricElement *ref
     ) const; 
    
  };
  
}

#endif

