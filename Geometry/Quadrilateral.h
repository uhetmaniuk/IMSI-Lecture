#ifndef FECODE_QUADRILATERAL_H
#define FECODE_QUADRILATERAL_H

#include <vector>

#include "Geometry/GeometricElement.h"

namespace FECode {
  
  class Vertex;
  
  //! A concrete class to define a quadrilateral in a mesh.
  
  /// The Quadrilateral class defines a quadrilateral with straight edges
  /// in a mesh on a domain.
  /// The class is derived from the virtual class GeometricElement.
  /// The vertices are numbered as follows
  ///
  ///    4 --- 3
  ///    |     |
  ///    1 --- 2
  ///
  class Quadrilateral: public virtual GeometricElement {
    
  protected:
    
  public:
    
    //////////////////////
    // ACCESSOR
    //////////////////////


    /// Returns the list of coordinates for all the vertices.
    std::vector<double> GetCoordinatesList() const;

    
    /// Computes the Jacobian of the map from the reference element
    /// to the physical element at a local point.
    /// The reference element is here the unit square [0,1]x[0,1].
    /// \param[in] lCoord: Vector of local coordinates of the point 
    /// in [0,1] x [0,1].
    /// \param[in,out] B: Jacobian matrix of the map at the local point.
    /// \return Boolean flag, set to TRUE when the computation is possible
    /// else set to FALSE.
    /// \note In 2D, the Jacobian matrix is [dx/di, dx/eta; dy/dxi, dy/deta].
    /// \note In 3D, the Jacobian matrix is not implemented (05/2008).
    bool GetMapJacobian
    (
     const std::vector<double> &lCoord,
     std::vector<double> &B
     ) const;
    
    /// Returns the measure of the element.
    double GetMeasure() const;
    
    /// Returns the list of side elements.
    /// \note In 2D, each edge [N_1, N_2] is such that rotating the vector N_1 N_2 
    /// by -90 degrees gives a normal outer to the triangular element.
    std::vector<GeometricElement*> GetSideList() const;
    
    /// Returns the number of vertices in a quadrilateral.
    int NumberVertices() const { return 4; };
    
    /// Returns the number of sides in a quadrilateral.
    int NumberSides() const { return 4; };
    
    ////////////
    // CREATOR
    ////////////
    
    /// Constructor.
    /// \param[in] dim - Space dimension (short).
    /// \param[in] vertexList - List of vertices forming the quadrilateral.
    /// \param[in] boundFlag - Boolean to flag whether the quadrilateral is on the boundary (Default value = 'false').
    Quadrilateral
    (
     short dim, 
     std::vector<Vertex*> &vertexList,
     bool boundFlag = false
     );
    
    /// Destructor.
    virtual ~Quadrilateral();
    
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

