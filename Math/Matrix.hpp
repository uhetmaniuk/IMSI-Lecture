#ifndef MATHCXX_MATRIX_HPP
#define MATHCXX_MATRIX_HPP


#include "Morphism.hpp"


namespace MathCXX 
{

  //! A templated abstract class to store a matrix.
  
  /// The Matrix class enables the use of matrices.
  template <class Scalar>
  class Matrix: public virtual Morphism<Scalar>
  {
 
    public:
   
 
      /////////////////
      // CREATOR
      ////////////////


      //! Destructor.
      virtual ~Matrix() { };


      ///////////////
      // MANIPULATOR
      ///////////////


      //! Add a value to an entry.
      //! \param[in] gRow  Index of the row where to add a scalar using C-notation.
      //! \param[in] gCol  Index of the column where to add a scalar using C-notation.
      //! \param[in] val  Scalar to add in the specified entry.
      //! \return Flag whether the operation was successful or not.
      virtual int Add
      (
        int gRow,
        int gCol,
        Scalar val
      )
      { return -1; };

    
      //! Add a block of values.
      //! \param[in] pos  List of row and column indices using C-notation.
      //! \param[in] squareMat  Submatrix stored in row-wise.
      //! \return Flag whether the operation was successful or not.
      virtual int Add
      (
        const std::vector<int> &pos,
        const std::vector<Scalar> &squareMat
      ) 
      { return -1; };


      /////////////////
      // ACCESSOR
      ////////////////


      /// Returns the number of rows.
      /// \return Number of rows.
      virtual size_t NumRows
      (
      ) const = 0;


      /// Returns the number of columns.
      /// \return Number of columns
      virtual size_t NumColumns
      (
      ) const = 0;

  };

}

#endif

