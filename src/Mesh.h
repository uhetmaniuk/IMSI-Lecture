#pragma once

#include "Element.h"

#include <array>
#include <ostream>
#include <vector>

namespace IMSI {
 
  /// \brief A class to manage the mesh on a domain.
  ///
  /// The Mesh class defines the data structure to manage a mesh on a bounded domain.
  class Mesh {
    
  protected:

    std::vector< std::array<double, 3> > vertexList;
    std::vector<ElementType> cellType;
    std::vector< std::vector<int> > cellToNode;

    /// \brief Array of node indices lying on the domain surface
    /// \note This structure should be generated for general boundary conditions.
    std::vector< int > boundaryNode;

  public:
    
    
    //////////////////////
    // ACCESSOR
    //////////////////////
    
    /// Returns the coordinates of a vertex.
    /// \param[in] id: ID number of the vertex.
    /// \return Array of coordinates
    std::array<double, 3> GetVertex
    (
     const int id
     ) const {
        return vertexList[id];
    }
    
    /// Returns the number of cells in the mesh.
    /// \return Number of cells.
    auto NumberCells() const 
    { return cellToNode.size(); };
    

    /// Returns the number of vertices in the mesh.
    /// \return Number of vertices.
    auto NumberVertices () const 
    { return vertexList.size(); };

    /// Returns the type of element in the specified cell
    /// \return Element type
    [[nodiscard]] ElementType GetCellType(int id) const {
        return cellType[id];
    }

    [[nodiscard]] const std::vector<std::vector<int>>& CellToNode() const {
        return cellToNode;
    }

    ////////////
    // CREATOR
    ////////////
    
    
    /// Constructor.
    /// \param[in] dim: Value of the spatial dimension (short integer).
    /// \param[in] nVertices: Number of vertices in the mesh (integer).
    /// \param[in] nCells: Number of cells in the mesh (integer).
    /// \note 'nVertices' and 'nCells' are optional parameters.
    /// When known, they allow to pre-allocate the memory for storing the vertices and cells.
    Mesh
    (
            std::vector< std::array<double, 3> >&& vertex,
      std::vector<ElementType>&& cList,
      std::vector< std::vector<int> >&& cToN,
      std::vector<int>&& nodeBdry
     ) : vertexList(std::move(vertex)), cellType(std::move(cList)), cellToNode(std::move(cToN)),
     boundaryNode(std::move(nodeBdry)) {}

    ////////////////
    // MANIPULATOR
    ////////////////
    
    //////////
    // FRIEND
    //////////
    
    
    /// Overloads the output operator <<.
    friend std::ostream & operator<<
    (
     std::ostream & os, 
     const IMSI::Mesh &ref
     );

  };
  
}

