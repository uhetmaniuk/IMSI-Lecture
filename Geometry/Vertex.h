#ifndef FECODE_VERTEX_H
#define FECODE_VERTEX_H

#include <vector>

#include "Geometry/GeometricElement.h"

namespace FECode {
  
  //! A concrete class to define a vertex in a mesh.
  
  /// The Vertex class defines a vertex in a mesh on a domain.
  /// The class is derived from the virtual class GeometricElement.
  /// The coordinates of the vertex are in the cartesian system.
  class Vertex: public virtual GeometricElement {
    
  private:
    
    // Don't define this function
    Vertex & operator=(const Vertex &ref);
    
  protected:
    
    std::vector<double> d_coordinates;
    
  public:
    
    //////////////////////
    // ACCESSOR
    //////////////////////
    
    /// Returns the list of coordinates.
    std::vector<double> GetCoordinates() const;
    
    /// Returns the list of coordinates.
    std::vector<double> GetCoordinatesList() const
    { return GetCoordinates(); };
    
    /// Returns the measure of the element.
    double GetMeasure() const { return 0.0; };
    
    /// Returns the list of side elements.
    std::vector<GeometricElement*> GetSideList() const;
    
    /// Returns the list of vertices.
    std::vector<Vertex*> GetVertexList() const;
    
    /// Returns the number of vertexs.
    int NumberVertices() const { return 1; };
    
    /// Returns the number of sides for a vertex.
    int NumberSides() const { return 0; };
    
    ////////////
    // CREATOR
    ////////////
    
    /// Constructor.
    /// \param[in] dim - Space dimension (short).
    /// \param[in] cood - List of coordinates for the vertex.
    /// \param[in] boundFlag - Boolean to flag whether the vertex is on the boundary (Default value = 'false').
    Vertex
    (
     short dim, 
     std::vector<double> &coord, 
     bool boundFlag = false
     );
    
    /// Constructor
    //!
    //! This function creates a hard copy of the reference vertex.
    Vertex(const Vertex &ref);
    
    /// Destructor.
    ~Vertex();
    
    ///////////////
    // MANIPULATOR
    ///////////////
    
    /// Checks whether the two elements are matching.
    /// \param[in] ref: Pointer to the element that neeeds to be compared.
    /// \return Boolean value whether the two elements match.
    bool CheckMatch(GeometricElement *ref) const; 
    
    /// Sets the coordinates.
    /// \param[in] coord: Vector of new coordinates for the vertex.
    int SetCoordinates(const std::vector<double> &coord);
    
    ///////////
    // FRIEND
    ///////////
    
    friend bool operator == (const FECode::Vertex &n1, const FECode::Vertex &n2);
    
  };
  
}

#endif

