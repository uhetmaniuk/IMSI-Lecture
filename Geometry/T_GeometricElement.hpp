#ifndef FECODE_T_GEOMETRIC_ELEMENT_HPP
#define FECODE_T_GEOMETRIC_ELEMENT_HPP


#include <iostream>
#include <vector>


namespace FECode 
{

  
  enum ElementType {VERTEX, EDGE, HEXAHEDRON, QUADRILATERAL, TETRAHEDRON, TRIANGLE};
  enum CoordinateType {CARTESIAN, CYLINDRICAL, POLAR, SPHERICAL};

  class Vertex;

  
  //! An abstract class to define a geometric element in a mesh.

  /// The T_GeometricElement class is a templated virtual class
  /// (specifies interface) that defines a geometric element
  /// in a mesh on a domain.
  template <class Scalar>
  class T_GeometricElement 
  {

    protected:

      short d_dim;

      bool d_onBoundary;

      ElementType d_Type;

      CoordinateType d_Coordinate;

      std::vector<Vertex*> d_vertexList;

      ///////////////////////
      // CREATOR (PROTECTED)
      ///////////////////////

      // Constructor
      T_GeometricElement<Scalar>
      (
        short dim
        , bool boundFlag
        , ElementType eType
        , CoordinateType cType = CARTESIAN
      );

    public:

      /////////////
      // ACCESSOR 
      /////////////

      /// Gets the boundary flag.
      bool GetBoundaryFlag() const
           { return d_onBoundary; };

      /// Gets the type of coordinate system.
      CoordinateType GetCoordinateType() const
           { return d_Coordinate; };

    /// Gets the type of coordinate system.
    virtual std::vector<Scalar> GetCoordinatesList() const = 0;
    
    /// Returns the spatial dimension.
      short GetDimension() const 
            { return d_dim; };

      /// Computes the Jacobian of the map from the reference element
      /// to the physical element at a local point.
      /// \param[in] lCoord Vector of local coordinates of the point. 
      /// \param[in,out] B Jacobian matrix of the map at the local point.
      /// \return Boolean flag, set to TRUE when the computation is possible
      /// else set to FALSE.
      virtual bool GetMapJacobian
      (
        const std::vector<Scalar> &lCoord,
        std::vector<Scalar> &B
      ) const;

      /// Returns the measure of the element (2D: area, 3D: volume).
      virtual Scalar GetMeasure() const = 0;

      /// Returns the list of vertices.
      virtual std::vector<Vertex*> GetVertexList() const;

      /// Returns the list of side elements.
      virtual std::vector< T_GeometricElement<Scalar>* > GetSideList() const = 0;

      /// Returns the type of the element.
      ElementType GetType() const 
            { return d_Type; };

      /// Returns the number of vertices for this element.
      virtual int NumberVertices() const = 0;

      /// Returns the number of sides for this element.
      virtual int NumberSides() const = 0;

      ////////////
      // CREATOR
      ////////////

      /// Virtual destructor.
      virtual ~T_GeometricElement() { };

      ///////////////
      // MANIPULATOR
      ///////////////

      /// Checks whether the two elements are matching.
      /// \param[in] ref Pointer to the element that neeeds to be compared.
      /// \return Boolean value whether the two elements match.
      virtual bool CheckMatch
      (
        T_GeometricElement<Scalar> *ref
      ) const = 0; 

      /// Sets the boundary flag.
      /// \param[in] bFlag Boolean determining whether the element
      /// is on the boundary.
      void SetBoundaryFlag(bool bFlag) 
           { d_onBoundary = bFlag; };

      /// Sets the type of coordinate system.
      /// \param[in] cFlag CoordinateType flag for the coordinate system.
      void SetCoordinateType(CoordinateType cType) 
           { d_Coordinate = cType; };

      /// Sets the spatial dimension.
      /// \param[in] dim Short integer to define the spatial dimension.
      void SetDimension(short dim) 
           { d_dim = dim; };

  };

  
  /////////////////////////////////////
  //
  // Definition of templated member functions
  //
  ////////////////////////////////////


  template<class Scalar>
  T_GeometricElement<Scalar>::T_GeometricElement
  (
   short dim
   , bool boundFlag
   , ElementType eType
   , CoordinateType cType
  ) :
    d_dim(dim)
  , d_onBoundary(boundFlag)
  , d_Type(eType)
  , d_Coordinate(cType)
  , d_vertexList()
  {
  }
  
  
  template<class Scalar>
  bool T_GeometricElement<Scalar>::GetMapJacobian
  (
   const std::vector<Scalar> &lCoord,
   std::vector<Scalar> &B
  ) const 
  {
    
    std::cerr << " !!! The computation of the Jacobian is not implemented !!! ";
    std::cerr << std::endl;
    
    return false;
    
  }


  template<class Scalar>
  std::vector<Vertex*> T_GeometricElement<Scalar>::GetVertexList
  (
  ) const
  {

    std::vector<Vertex*> list(d_vertexList.begin(), d_vertexList.end());
    return list;

  }


}

#endif

