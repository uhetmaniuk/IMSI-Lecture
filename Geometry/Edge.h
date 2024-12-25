#ifndef FECODE_EDGE_H
#define FECODE_EDGE_H

#include <vector>

#include "Geometry/GeometricElement.h"

namespace FECode {

  class Vertex;

  //! A concrete class to define an edge in a mesh.

  /// The Edge class defines an edge in a mesh on a domain.
  /// The class is derived from the virtual class GeometricElement.
  class Edge: public virtual GeometricElement {

    protected:

    public:

      //////////////////////
      // ACCESSOR
      //////////////////////

      /// Returns the list of coordinates.
      std::vector<double> GetCoordinatesList() const;
    
      /// Returns the measure of the element.
      double GetMeasure() const;
 
      /// Returns the list of side elements.
      std::vector<GeometricElement*> GetSideList() const;

      /// Returns the number of vertices in an edge.
      int NumberVertices() const { return 2; };

      /// Returns the number of sides in an edge.
      int NumberSides() const { return 2; };

      ////////////
      // CREATOR
      ////////////

      /// Constructor.
      /// \param[in] dim - Space dimension (short).
      /// \param[in] vertexList - List of vertices forming the edge.
      /// \param[in] boundFlag - Boolean to flag whether the edge is on the boundary (Default value = 'false').
      Edge
      (
        short dim, 
        std::vector<Vertex*> &vertexList,
        bool boundFlag = false
      );

      /// Destructor.
      virtual ~Edge();

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

