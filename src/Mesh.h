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

      int sdim;

      std::vector< double > vertex_x;
      std::vector< double > vertex_y;
      std::vector< double > vertex_z;
    std::vector<ElementType> cellType;
    std::vector< std::vector<int> > cellToNode;

    /// \brief Array of node indices lying on the domain surface
    /// \note This structure should be generalized for other boundary conditions.
    std::vector< int > boundaryNode;

  public:


      /// Constructor.
      Mesh(
                      int dim,
                      std::vector< std::array<double, 3> > const& vertex,
                      std::vector<ElementType>&& cList,
                      std::vector< std::vector<int> >&& cToN,
                      std::vector<int>&& nodeBdry
              );

    /// Returns the coordinates of a vertex.
    /// \param[in] id: ID number of the vertex.
    /// \return Array of coordinates
    [[nodiscard]] std::array<double, 3> GetVertex
    (
     const int id
     ) const {
        return {vertex_x[id], vertex_y[id], vertex_z[id]};
    }
    
    /// Returns the number of cells in the mesh.
    /// \return Number of cells.
    [[nodiscard]] auto NumberCells() const
    { return cellToNode.size(); };
    

    /// Returns the number of vertices in the mesh.
    /// \return Number of vertices.
    [[nodiscard]] auto NumberVertices () const
    { return vertex_x.size(); };

    /// Returns the type of element in the specified cell
    /// \return Element type
    [[nodiscard]] ElementType GetCellType(int id) const {
        return cellType[id];
    }

    [[nodiscard]] const std::vector<std::vector<int>>& CellToNode() const {
        return cellToNode;
    }

      [[nodiscard]] const std::vector<int>& NodeList(int cellID) const {
          return cellToNode[cellID];
      }

      [[nodiscard]] auto GetSpatialDimension() const {
        return sdim;
    }

    [[nodiscard]] const std::vector<int>& GetBoundaryNodes() const {
        return boundaryNode;
    }

    /// Overloads the output operator <<.
    friend std::ostream & operator<<
    (
     std::ostream & os, 
     const IMSI::Mesh &ref
     );

  };
  
}

