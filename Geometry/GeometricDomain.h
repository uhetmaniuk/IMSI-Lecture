#ifndef FECODE_GEOMETRIC_DOMAIN_H
#define FECODE_GEOMETRIC_DOMAIN_H


#include <string>
#include <vector>


#ifdef FEC_USE_TRILINOS
class Epetra_Map;
#endif


namespace MathCXX {

  class Graph;

}


namespace FECode {

  
  class Communicator;
  class Mesh;
  

  //! An abstract class to define a geometric domain.
  
  /// The GeometricDomain class is a virtual class
  /// (specifies interface) that defines a geometric domain.
  class GeometricDomain {
    
  private:
    
    // Do not define these functions.
    
    GeometricDomain(const GeometricDomain &ref);
    GeometricDomain & operator=(const GeometricDomain &ref);
    
  protected:
    
    // Define the data.
    const Communicator &d_Comm;
    
    short d_dim;

    Mesh *d_Mesh_p;
    
    // (UH) 05/10
    // The pointer d_GlobalVertexToVertex_p is used to keep track of the global mesh.
    // It works for nobal-based finite elements.
    // For edge-based finite elements, we might need something else (cell-to-cell ?)
    // We might consider using Epetra_CrsGraph.
    MathCXX::Graph *d_Global_CellToCell_p;
    MathCXX::Graph *d_Global_VertexToVertex_p;
    
#ifdef FEC_USE_TRILINOS
    Epetra_Map *d_CellMap_p;
    Epetra_Map *d_VertexMap_p;
    Epetra_Map *d_VertexMap_1to1_p;
#endif
    
    std::string d_domainLabel;
    
    
    ///////////////////////
    // CREATOR (PROTECTED)
    ///////////////////////
    
    
    // Constructor
    GeometricDomain
    (
     const Communicator &Comm,
     short dim
     );
    
    
    ///////////////////////////
    // MANIPULATOR (PROTECTED)
    ///////////////////////////
    
    
    // Partitions the domain into subdomains.
    // The default algorithms use METIS or a uniform distribution.
    void DefaultPartitionDomain
    (
      std::vector<int> &cellToPart
    ) const;
    
    
    // Distributes the global mesh into one local mesh per subdomain.
    // Note: the global mesh will be deleted.
    void DistributeGlobalMesh
    (
      const std::vector<int> &cellToPart
     );
    
    
    // Partitions the domain into subdomains.
    virtual void PartitionDomain();
    
    
  public:
    

    /////////////
    // ACCESSOR 
    /////////////
    
    
#ifdef FEC_USE_TRILINOS
    /// Returns the one-to-one map of vertices among processors.
    /// \return Pointer to Epetra_Map object.
    Epetra_Map* Get_1to1_VertexMap
    (
    ) const
    {
      return d_VertexMap_1to1_p;
    };
#endif

 
    /// Returns the pointer to the mesh.
    /// \return Pointer to Mesh object.
    Mesh* GetMesh
    (
    ) const
    {
      return d_Mesh_p;
    };
    

    /// Returns the pointer to the connectivity between vertices.
    /// \return Pointer to Graph object.
    MathCXX::Graph* Get_Global_VertexToVertex
    (
    ) const;
    
    
    /// Returns a label describing the domain.
    virtual const std::string Label
    (
    ) const;

    
    /// Returns the global number of vertices.
    /// \return Number of vertices.
    /// \note The 'int' might get too small for large meshes.
    int NumberGlobalVertices
    (
    ) const;

    
    /// Returns the spatial dimension of the mesh.
    /// \return Spatial dimension.
    short int SpatialDim
    (
    ) const 
    { 
    return d_dim; 
    };
    

    ////////////
    // CREATOR
    ////////////
    
    
    /// Virtual destructor.
    virtual ~GeometricDomain();
    
    
    ///////////////
    // MANIPULATOR
    ///////////////
    
    
  };
  
}

#endif

