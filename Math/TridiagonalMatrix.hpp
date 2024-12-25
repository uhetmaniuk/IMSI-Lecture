#ifndef MATH_CXX_TRIDIAGONAL_MATRIX_HPP
#define MATH_CXX_TRIDIAGONAL_MATRIX_HPP 


#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <sys/time.h>
#include <vector>


#include "BlockVector.hpp"
#include "Matrix.hpp"


namespace MathCXX {

  
  //! A templated class for tridiagonal matrices.
  
  /// The TridiagonalMatrix class enables the use of tridiagonal matrices.
  /// \note For a nxn matrix, the first n values correspond to diagonal entries.
  /// The next (n-1) values correspond to the upper diagonal entries.
  /// The next (n-1) values are for the lower diagonal entries.
  /// \note (2011/01) The implementation is only for square matrices.
  template <class Scalar>
  class TridiagonalMatrix: public virtual Matrix<Scalar> {
    
  private:
    
    
  protected:


    size_t d_numRow;
    size_t d_numCol;
    Scalar *d_values_p;


  public:


    ///////////
    // CREATOR
    ///////////
    
    
    //! Constructor.
    //! \param[in] nRow  Number of rows.
    //! \param[in] nCol  Number of columns.
    /// \note The values of the matrix are not initialized.
    /// \note (2011/01) The implementation is only for square matrices.
    TridiagonalMatrix<Scalar>
    (
      size_t nRow,
      size_t nCol
     );
    
    
    //! Copy constructor.
    //! \param[in] ref Matrix to copy.
    /// \note Only the dimensions and the values are copied.
    TridiagonalMatrix<Scalar>
    (
     const TridiagonalMatrix<Scalar> &ref
     );
    
    
    //! Destructor.
    virtual ~TridiagonalMatrix
    (
    );
    
    
    ///////////////
    // MANIPULATOR 
    ///////////////
    
    
    //! Set all values in a matrix with constant value.
    /// \param[in] arg Value to use.
    /// \return Updated matrix-vector.
    TridiagonalMatrix<Scalar>& operator=
    (
     const Scalar &arg
     );


    //! Add a value to an entry.
    //! \param[in] gRow Index of the row where to add a scalar using C-notation.
    //! \param[in] gCol Index of the column where to add a scalar using C-notation.
    //! \param[in] val Scalar to add in the specified entry.
    //! \return Flag whether the operation was successful or not.
    int Add
    (
      int gRow,
      int gCol,
      Scalar val
    );

    
    //! Add a block of values.
    //! \param[in] pos List of row and column indices using C-notation.
    //! \param[in] squareMat Submatrix stored in row-wise.
    //! \return Flag whether the operation was successful or not.
    //! \note The routine does not do any check about squareMat.
    //! It adds only the values for the tridiagonal of the matrix.
    int Add
    (
      const std::vector<int> &pos,
      const std::vector<Scalar> &squareMat
    ); 


    //! Apply the matrix to a block of vectors
    //! \param[in] X Block vectors on which the matrix is applied.
    //! \param[in] Y Block vectors storing the result.
    int Apply
    (
      const BlockVector<Scalar> &X,
      BlockVector<Scalar> &Y
    );


    //! Sets multi-vector values to random numbers.
    /// \param[in] seed  Integer setting the state of the random generator.
    /// \note A seed negative or zero will reset the state to the time.
    /// \note The real entries will be between -1 and 1.
    void Random
    (
     int seed = 0
    );
    
    
    ////////////
    // ACCESSOR 
    ////////////
    
   
    //! Returns the number of rows in the matrix.
    size_t NumRows
    () const
    { return d_numRow; };
    
    
    //! Returns the number of columns in the matrix.
    size_t NumColumns
    () const
    { return d_numCol; };
    
    
    //! Returns a reference to the pointer of values.
    /// \note This function gives a read-only access to the pointer.
    /// \note For a nxn matrix, the first n values correspond to diagonal entries.
    /// The next (n-1) values correspond to the upper diagonal entries.
    /// The next (n-1) values are for the lower diagonal entries.
    Scalar* Values
    () const
    { 
      return d_values_p;
    };
    
    
    /// Element access function
    /// \param[in] RowIndex Index of row for the element using C-notation.
    /// \param[in] ColIndex Index of column for the element using C-notation.
    /// \return Element from the specified row and column.
    /// \warning No bounds checking is done.
    /// \note This function allows write access to a value.
    Scalar& operator() 
    (
     size_t RowIndex, 
     size_t ColIndex
     );
    
    
    /// Element access function
    /// \param[in] RowIndex Index of row for the element using C-notation.
    /// \param[in] ColIndex Index of column for the element using C-notation.
    /// \return Element from the specified row and column.
    /// \warning No bounds checking is done.
    /// \note This function is similar to a read-only access.
    Scalar operator() 
    (
     size_t RowIndex, 
     size_t ColIndex
     ) 
    const;
    
    
    //////////
    // FRIEND
    //////////
    
    
  };
  
  
}


//==============================================================================
//==============================================================================


namespace MathCXX {
  
  
  /////////////////////////////////////
  //
  // Definition of templated member functions
  //
  ////////////////////////////////////
  
  
  template<class Scalar>
  TridiagonalMatrix<Scalar>::TridiagonalMatrix
  (
   size_t nRow,
   size_t nCol
   )
  :
  d_numRow(nRow),
  d_numCol(nCol),
  d_values_p(new Scalar[3*nRow-1])
  {
    // The array d_values_p has one optional entry.
    // It is allocated to avoid returning a local variable with (.,.)
    assert(nCol == nRow);
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  TridiagonalMatrix<Scalar>::TridiagonalMatrix
  (
   const TridiagonalMatrix<Scalar> &ref
   )
  :
  d_numRow(ref.d_numRow),
  d_numCol(ref.d_numCol),
  d_values_p(new Scalar[3*d_numRow-1])
  {
    
    // The array d_values_p has one optional entry.
    // It is allocated to avoid returning a local variable with (.,.)

    size_t len = 3*d_numRow - 1; 
    for (size_t ii = 0; ii < len; ++ii)
      d_values_p[ii] = ref.d_values_p[ii];
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  TridiagonalMatrix<Scalar>::~TridiagonalMatrix
  ()
  {
    if (d_values_p)
      delete[] d_values_p;
    d_values_p = 0;
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  TridiagonalMatrix<Scalar>& TridiagonalMatrix<Scalar>::operator=
  (
   const Scalar &arg
   )
  {
    size_t len = 3*d_numRow - 1;
    for (size_t ii = 0; ii < len; ++ii)
      d_values_p[ii] = arg;
    return *this;
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  int TridiagonalMatrix<Scalar>::Add
  (
   int gRow,
   int gCol,
   Scalar val
   )
  {
    
    if (gRow == gCol)
    {
      // Diagonal entry
      d_values_p[gRow] += val;
      return 0;
    }
    else if (gRow == gCol - 1)
    {
      // Upper diagonal entry
      d_values_p[d_numRow + gRow] += val;
      return 0;
    }
    else if (gRow == gCol + 1)
    {
      // Lower diagonal entry
      d_values_p[2*d_numRow + gCol - 1] += val;
      return 0;
    }
    
    return -1;
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  int TridiagonalMatrix<Scalar>::Add
  (
   const std::vector<int> &pos,
   const std::vector<Scalar> &squareMat
   )
  {
    
    int info = -1;
    
    bool isGood = true;
    
    size_t iLen = pos.size();
    if (squareMat.size() != iLen * iLen)
      return -1;
    
    for (size_t ii = 0; ii < iLen; ++ii)
    {
      if ((pos[ii] < 0) || (pos[ii] >= d_numRow))
      {
        isGood = false;
        break;
      }
    }
    
    if (isGood == false)
      return info;
    
    for (size_t ii = 0; ii < iLen; ++ii)
    {
      int iRow = pos[ii];
      for (size_t jj = 0; jj < iLen; ++jj)
      {
        int jCol = pos[jj];
        if (iRow == jCol)
        {
          // Diagonal entry
          d_values_p[iRow] += squareMat[ii + jj * iLen];
        }
        else if (iRow == jCol - 1)
        {
          // Upper diagonal entry
          d_values_p[d_numRow + iRow] += squareMat[ii + jj * iLen];
        }
        else if (iRow == jCol + 1)
        {
          // Lower diagonal entry
          d_values_p[2*d_numRow - 1 + jCol] += squareMat[ii + jj * iLen];
        }
      } // for (size_t jj = 0; jj < iLen; ++jj)
    }
    
    return info;
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  int TridiagonalMatrix<Scalar>::Apply
  (
   const BlockVector<Scalar> &X
   , BlockVector<Scalar> &Y
   )
  {
    
    int info = -1;
    
    if ((X.NumRows() != d_numCol) || (Y.NumRows() != d_numRow))
      return -1;
    
    size_t numVec = X.NumVectors();
    if (numVec > Y.NumVectors())
      return -1;
    
    Scalar *Xval_p = X.Values();
    Scalar *Yval_p = Y.Values();
    
    for (size_t ii = 0; ii < d_numRow; ++ii)
    {
      
      Scalar coeff = d_values_p[ii];
      for (size_t jj = 0; jj < numVec; ++jj)
        Yval_p[ii + jj * d_numRow] = coeff * Xval_p[ii + jj * d_numCol];
      
      if (ii + 1 < d_numRow)
      {
        coeff = d_values_p[d_numRow + ii];
        for (size_t jj = 0; jj < numVec; ++jj)
          Yval_p[ii + jj * d_numRow] += coeff * Xval_p[ii+1 + jj * d_numCol];
      }

      if (ii > 0)
      {
        coeff = d_values_p[2*d_numRow - 1 + ii - 1];
        for (size_t jj = 0; jj < numVec; ++jj)
          Yval_p[ii + jj * d_numRow] += coeff * Xval_p[ii-1 + jj * d_numCol];
      }

    }
    
    return info;
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  void TridiagonalMatrix<Scalar>::Random
  (
   int seed
   )
  {
    
    size_t len = 3 * d_numRow - 1;
    for (size_t ii = 0; ii < len; ++ii)
      SetToRandom(d_values_p[ii]);
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  Scalar& TridiagonalMatrix<Scalar>::operator()
  (
   size_t RowIndex
   , size_t ColIndex
   )
  {
    
    if (RowIndex == ColIndex)
    {
      // Diagonal entry
      return d_values_p[RowIndex];
    }
    else if (RowIndex == ColIndex - 1)
    {
      // Upper diagonal entry
      return d_values_p[d_numRow + RowIndex];
    }
    else if (RowIndex == ColIndex + 1)
    {
      // Lower diagonal entry
      return d_values_p[2*d_numRow - 1 + ColIndex];
    }

    d_values_p[3*d_numRow - 2] = 0;
    return d_values_p[3*d_numRow - 2];
    
  }
  
  
  //-------------------------------------------------
  
  
  template<class Scalar>
  Scalar TridiagonalMatrix<Scalar>::operator()
  (
   size_t RowIndex
   , size_t ColIndex
   ) const
  {
    
    Scalar tmpval = 0;
    
    if (RowIndex == ColIndex)
    {
      // Diagonal entry
      tmpval = d_values_p[RowIndex];
    }
    else if (RowIndex == ColIndex - 1)
    {
      // Upper diagonal entry
      tmpval = d_values_p[d_numRow + RowIndex];
    }
    else if (RowIndex == ColIndex + 1)
    {
      // Lower diagonal entry
      tmpval = d_values_p[2*d_numRow - 1 + ColIndex];
    }
    
    return tmpval;
    
  }
  
  
  //-------------------------------------------------
  
  
}


#endif
