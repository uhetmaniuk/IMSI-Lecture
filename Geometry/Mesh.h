#pragma once

#include <vector>
#include <ostream>

#include "Geometry/GeometricElement.h"

namespace msfem {

class Graph;
  
  class Vertex;
  
  //! A class to manage the mesh on a domain.
  
  /// The Mesh class defines the data structure to manage a mesh on a bounded domain.
  class Mesh {
    
    
  protected:
        
    std::vector<real_t> d_nodeData;
      
    std::vector<GeometricElement*> d_cellList;
  
    public:
      
      using real_t = float;
    
      static constexpr char d_dim = 3;

      Mesh() = delete;
      
      Mesh(size_t nNodes, size_t nCells) : d_nodeData(d_dim * nNodes, 0.0),
      d_cellList(nCells, nullptr) {}
      
      ~Mesh();
      
      Mesh& operator=(const Mesh &ref) = delete;
      
      size_t Nodes() const { return d_nodeList.size() / d_dim; }
      
      size_t Cells() const { return d_cellList.size(); }
      
      void SetNode(size_t id, Span<real_t> coord);
      

      /*
    //////////////////////
    // ACCESSOR
    //////////////////////
    
    
    /// Returns a pointer to a cell.
    /// \param[in] id: ID number of the cell.
    /// \return Pointer to the cell, set to 0 if the search is not successful.
    GeometricElement* GetCell
    (
     const int id
     ) const;
    
    
    /// Returns the ID number for a cell.
    /// \param[in] thisCell: Pointer to the geometric element.
    /// \return Integer ID number, set to -1 if the cell is not stored.
    int GetCellID
    (
     const GeometricElement *thisCell_p
     ) const;
    
    
    /// Returns a pointer to a vertex.
    /// \param[in] id: ID number of the vertex.
    /// \return Pointer to the vertex, set to 0 if the search is not successful.
    Vertex* GetVertex
    (
     const int id
     ) const;
    
    
    /// Returns the ID number for a vertex.
    /// \param[in] thisVertex_p: Pointer to the vertex.
    /// \return Integer ID number, set to -1 if the vertex is not stored.
    int GetVertexID
    (
     const Vertex *thisVertex_p
     ) const;
    
    
    /// Returns a pointer to a side.
    /// \param[in] id: ID number of the side.
    /// \return Pointer to the side, set to 0 if the search is not successful.
    GeometricElement* GetSide
    (
     const int id
     ) const;
    
    
    /// Returns the number of cells in the mesh.
    /// \return Number of cells.
    int NumberCells
    (
    ) const 
    { return d_numCells; };
    

    /// Returns the number of vertices in the mesh.
    /// \return Number of vertices.
    int NumberVertices
    (
    ) const 
    { return d_numVertices; };
   
 
    /// Returns the number of vertices in the mesh.
    /// \return Number of vertices.
    int NumVertices
    (
    ) const 
    { return d_numVertices; };
   
 
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
     int nVertices = 0, 
     int nCells = 0
     );
    
    
    /// Constructor
    //!
    //! This function creates a hard copy of the reference mesh.
    Mesh
    (
     const Mesh &ref
     );
    
    
    /// Destructor.
    virtual ~Mesh();
    
    
    ////////////////
    // MANIPULATOR
    ////////////////
    
    
    /// Inserts a new vertex to mesh.
    /// \param[in] newVertex_p: Pointer to the new vertex.
    /// \return Boolean result, set to TRUE if the insertion was successful 
    /// and FALSE if the pointer or the vertex is already stored.
    bool InsertVertex
    (
     Vertex *newVertex_p
     );
    
    
    /// Inserts a new cell to mesh.
    /// \param[in] newCell: Pointer to the new cell.
    /// \return Boolean result, set to TRUE if the insertion was successful 
    /// and FALSE if the pointer or the cell is already stored.
    bool InsertCell
    (
     GeometricElement *newCell_p
     );
    
    
    /// Makes the cell-to-cell connectivity in terms of indices.
    /// \param[in,out] cellToCell: Reference to a 'Graph' object.
    /// \param[in] useVertexConnectivity: Optional flag to use the nodal connectivity (default value = false).
    /// \note The cell-to-cell connectivity differs when using the nodal or side connectivity as
    /// as a step
    /// The default approach is to use the side connectivity.
    /// When 
    void MakeCellToCell
    (
     MathCXX::Graph &cellToCell,
     bool useVertexConnectivity = false
     ) const;
    
    
    /// Makes the cell-to-vertex connectivity in terms of indices.
    /// \param[in,out] cellToVertex: Reference to a 'Graph' object.
    void MakeCellToVertex
    (
     MathCXX::Graph &cellToVertex
     ) const;
    
    
    /// Makes the vertex-to-cell connectivity in terms of indices.
    /// \param[in,out] vertexToCell: Reference to a 'Graph' object.
    void MakeVertexToCell
    (
     MathCXX::Graph &vertexToCell
     ) const;
    
    
    /// Makes the vertex-to-vertex connectivity in terms of indices.
    /// \param[in,out] vertexToVertex: Reference to a 'Graph' object.
    void MakeVertexToVertex
    (
     MathCXX::Graph &vertexToVertex
     ) const;
    
    
    /// Makes the side-to-cell connectivity in terms of indices.
    /// \param[in,out] sideToCell: Reference to a 'Graph' object.
    /// \note This routine creates also the list of sides.
    void MakeSideToCell
    (
     MathCXX::Graph &sideToCell
     ) const;
    
    
    /// Set a specified cell.
    /// \param[in] cell_p: Pointer to the new cell.
    /// \param[in] id: ID number of the cell.
    /// \return Boolean  result, set to TRUE if the change was successful 
    /// and FALSE if not.
    bool SetCell
    (
     GeometricElement* cell_p, 
     int id
     );
    
    
    /// Set a specified side.
    /// \param[in] side_p: Pointer to the new side.
    /// \param[in] id: ID number of the side.
    /// \return Boolean  result, set to TRUE if the change was successful 
    /// and FALSE if not.
    bool SetSide
    (
     GeometricElement* side_p, 
     int id
     );
    
    
    /// Set a specified vertex.
    /// \param[in] cell_p: Pointer to the new vertex.
    /// \param[in] id: ID number of the vertex.
    /// \return Boolean  result, set to TRUE if the change was successful 
    /// and FALSE if not.
    bool SetVertex
    (
     Vertex* vertex_p, 
     int id
     );
    */
    
    //////////
    // FRIEND
    //////////
    
    
    /// Overloads the output operator <<.
    friend std::ostream & operator<<
    (
     std::ostream & os, 
     const msfem::Mesh &ref
     );
    
    
  };
  
} // namespace msfem
