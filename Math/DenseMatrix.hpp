#ifndef MATH_CXX_DENSEMATRIX_HPP
#define MATH_CXX_DENSEMATRIX_HPP 


#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sys/time.h>
#include <vector>


#include "BlockVector.hpp"
#include "Matrix.hpp"


namespace MathCXX {


  //! A templated class for dense matrices.
  
  /// The DenseMatrix class enables the use of dense matrices.
  template <class Scalar>
  class DenseMatrix: public virtual Matrix<Scalar> {
    
  private:
    
    
    template <class T>
    bool HasMatchingDimensions
    (
     const DenseMatrix<T> &xx
     ) const;
    
    
  protected:
    
    
    size_t d_numRow;
    size_t d_numCol;
    
    size_t d_leadDim;
    
    Scalar *d_val_p;
    bool d_hasAllocatedVal;
    
    Scalar *d_work_p;
    size_t d_lwork;
    
    
  public:
    
    
    ///////////
    // CREATOR
    ///////////
    
    
    //! Constructor.
    //! \param[in] nRow Number of rows.
    //! \param[in] nCol Number of columns.
    /// \note The values of the matrix are not initialized.
    DenseMatrix<Scalar>
    (
     size_t nRow, 
     size_t nCol
     );
    
    
    //! Constructor.
    //! \param[in] nRow Number of rows.
    //! \param[in] nCol Number of columns.
    //! \param[in] ldw Leading dimension for a column.
    /// \note The values of the matrix are not initialized.
    DenseMatrix<Scalar>
    (
     size_t nRow, 
     size_t nCol,
     size_t ldw
     );
    
    
    //! Constructor.
    //! \param[in] nRow Number of rows.
    //! \param[in] nCol Number of columns.
    //! \param[in] val_p Pointer for the data.
    /// \note The data are not copied.
    DenseMatrix<Scalar>
    (
     size_t nRow, 
     size_t nCol,
     Scalar *val_p
     );
    
    
    //! Constructor.
    //! \param[in] nRow Number of rows.
    //! \param[in] nCol Number of columns.
    //! \param[in] val_p Pointer for the data.
    //! \param[in] ldw Leading dimension for a row in val_p.
    /// \note The data are not copied.
    DenseMatrix<Scalar>
    (
     size_t nRow, 
     size_t nCol,
     size_t ldw,
     Scalar *val_p
     );
    
    
    //! Copy constructor.
    //! \param[in] ref Matrix to copy.
    /// \note Only the dimensions and the values are copied.
    DenseMatrix<Scalar>
    (
     const DenseMatrix<Scalar> &ref
     );
    
    
    //! Copy constructor.
    //! \param[in] ref Matrix to copy.
    /// \note The type 'T' for matrix  must allow conversion to 'Scalar' type.
    /// \note Only the dimensions and the converted values are copied.
    template<class T>
    DenseMatrix<Scalar>
    (
     const DenseMatrix<T> &ref
     );
    
    
    //! Destructor.
    virtual ~DenseMatrix
    (
    );
    
    
    ///////////////
    // MANIPULATOR 
    ///////////////
    
    
    //! Copy assignment.
    /// \param[in] bv Matrix used overwrite the dimensions and data.
    /// \return Updated matrix.
    DenseMatrix<Scalar>& operator=
    (
     const DenseMatrix<Scalar> &bv
     );
    
    
    //! Copy assignment.
    /// \param[in] bv Matrix used overwrite the dimensions and data.
    /// \return Updated matrix.
    /// \note The type 'T' for bv must allow conversion to the 'Scalar' type.
    template <class T>
    DenseMatrix<Scalar>& operator=
    (
     const DenseMatrix<T> &bv
     );
    
    
    //! Set all values in a matrix with constant value.
    /// \param[in] arg Value to use.
    /// \return Updated matrix-vector.
    DenseMatrix<Scalar>& operator=
    (
     const Scalar &arg
     );
    
    
    //! Componentwise multiplication with a scalar.
    /// \param[in] sb Scalar value used to multiply all the components.
    /// \return Updated block vectors.
    DenseMatrix<Scalar>& operator*=
    (
     const Scalar &b
     );
    
    
    //! Componentwise division by a scalar.
    /// \param[in] sb Scalar value used to divide all the components.
    /// \return Updated block vectors.
    DenseMatrix<Scalar>& operator/=
    (
     const Scalar &b
     );
    
    
    //! Componentwise substraction with a blockvector. 
    /// \param[in] bv Multivector used to do the substraction.
    /// \return Updated block vectors.
    /// \note The dimensions must match.
    /// \note The type 'T' for bv must allow conversion to the 'Scalar' type.
    template <class T>
    DenseMatrix<Scalar>& operator-=
    (
     const DenseMatrix<T> &bv
     );
    
    
    //! Componentwise substraction with a scalar.
    /// \param[in] sb Scalar value used to substract to all the components.
    /// \return Updated block vectors.
    DenseMatrix<Scalar>& operator-=
    (
     const Scalar &b
     );
    
    
    //! Componentwise addition with a blockvector. 
    /// \param[in] bv Multivector used to do the addition.
    /// \return Updated block vectors.
    /// \note The dimensions must match.
    /// \note The type 'T' for bv must allow conversion to the 'Scalar' type.
    template <class T>
    DenseMatrix<Scalar>& operator+=
    (
     const DenseMatrix<T> &bv
     );
    
    
    //! Componentwise addition with a scalar.
    /// \param[in] sb  Scalar value used to add to all the components.
    /// \return Updated block vectors.
    DenseMatrix<Scalar>& operator+=
    (
     const Scalar &b
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
    int Add
    (
      const std::vector<int> &pos,
      const std::vector<Scalar> &squareMat
    ); 


    //! Apply the matrix to a block of vectors
    //! \param[in] X Block vectors on which the matrix is applied.
    //! \param[in] Y Block vectors storing the result.
    virtual int Apply
    (
      const BlockVector<Scalar> &X,
      BlockVector<Scalar> &Y
    );
    
    
    //! Compute the weighted sum Y <- a*X + b*Y.
    /// \param[in] a  Scalar value for X.
    /// \param[in] X  Matrix to add.
    /// \param[in] b  Scalar value for Y.
    /// \return Updated matrix.
    /// \note The vector Y is given by (*this).
    /// \note The type 'T' for matrix X should allow conversion to 'Scalar' type for matrix Y.
    template <class T>
    void AXPBY
    (
     const Scalar &a,
     const DenseMatrix<T> &X,
     const Scalar &b
     );
    
    
    //! Compute the weighted sum Y(row, :) <- a*X + b*Y(row, :).
    /// \param[in] RowIndex  List of row indices.
    /// \param[in] a  Scalar value for X.
    /// \param[in] X  Matrix to add.
    /// \param[in] b  Scalar value for Y.
    /// \return Updated matrix.
    /// \note The vector Y is given by (*this).
    /// \note The type 'T' for matrix X should allow conversion to 'Scalar' type for matrix Y.
    /// \note The 'Index' type should allow conversion to 'size_t' type.
    template<class T, class Index>
    void local_AXPBY
    (
     const std::vector<Index> &RowIndex,
     const Scalar &a,
     const DenseMatrix<T> &X,
     const Scalar &b
     );
    
    
    //! Compute the weighted sum Y(row, col) <- a*X + b*Y(row, col).
    /// \param[in] RowIndex  List of row indices.
    /// \param[in] ColIndex  List of column indices.
    /// \param[in] a  Scalar value for X.
    /// \param[in] X  Matrix to add.
    /// \param[in] b  Scalar value for Y.
    /// \return Updated matrix.
    /// \note The multivector Y is given by (*this).
    /// \note The type 'T' for matrix X should allow conversion to 'Scalar' type for matrix Y.
    /// \note The 'Index' type should allow conversion to 'size_t' type.
    template<class T, class Index>
    void local_AXPBY
    (
     const std::vector<Index> &RowIndex,
     const std::vector<Index> &ColIndex,
     const Scalar &a,
     const DenseMatrix<T> &X,
     const Scalar &b
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
    
    
    //! Returns the leading dimension of the matrix.
    size_t LeadingDimension
    () const 
    { return d_leadDim; };
    
    
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
    Scalar* Values
    () const
    { return this->d_val_p; };
    
    
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
     )
    { return ( d_val_p[RowIndex + ColIndex * d_leadDim] ); };
    
    
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
    const
    { return ( d_val_p[RowIndex + ColIndex * d_leadDim] ); };
    
    
    /// Block access function
    /// \param[in] RowIndex Indices of rows using C-notation.
    /// \param[in] ColIndex Indices of columns using C-notation.
    /// \return Block from the specified rows and columns.
    /// \note This function is similar to a read-only access.
    /// \note The 'Index' type should allow conversion to 'size_t' type.
    template<class Index>
    const DenseMatrix<Scalar> operator() 
    (
     const std::vector<Index> &RowIndex, 
     const std::vector<Index> &ColIndex
     ) 
    const;
    
    
    //////////
    // FRIEND
    //////////
   

  };

  
}


//==============================================================================
//==============================================================================


#include "MathCXXCore.hpp"


namespace MathCXX {
  

  /////////////////////////////////////
  //
  // Definition of templated member functions
  //
  ////////////////////////////////////
  
  
  template<class Scalar>
  DenseMatrix<Scalar>::DenseMatrix
  (
   size_t nRow, 
   size_t nCol
   )
  : 
  d_numRow(nRow),
  d_numCol(nCol), 
  d_leadDim(d_numRow),
  d_val_p(),
  d_work_p(),
  d_lwork(0)
  {
    
    d_val_p = new Scalar[d_leadDim * d_numCol];
    d_hasAllocatedVal = true;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  DenseMatrix<Scalar>::DenseMatrix
  (
   size_t nRow, 
   size_t nCol,
   size_t ldw
   )
  : 
  d_numRow(nRow),
  d_numCol(nCol), 
  d_leadDim(ldw),
  d_val_p(),
  d_work_p(),
  d_lwork(0)
  {
    
    d_val_p = new Scalar[d_leadDim * d_numCol];
    d_hasAllocatedVal = true;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  DenseMatrix<Scalar>::DenseMatrix
  (
   size_t nRow,
   size_t nCol,
   Scalar *val_p
   )
  :
  d_numRow(nRow),
  d_numCol(nCol),
  d_leadDim(d_numRow),
  d_val_p(val_p),
  d_work_p(),
  d_lwork(0)
  {
    
    d_hasAllocatedVal = false;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  DenseMatrix<Scalar>::DenseMatrix
  (
   size_t nRow,
   size_t nCol,
   size_t ldw,
   Scalar *val_p
   )
  :
  d_numRow(nRow),
  d_numCol(nCol),
  d_leadDim(ldw),
  d_val_p(val_p),
  d_work_p(),
  d_lwork(0)
  {
    
    d_hasAllocatedVal = false;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>::DenseMatrix
  (
   const DenseMatrix<Scalar> &ref
   )
  :
  d_numRow(ref.NumRows()),
  d_numCol(ref.NumColumns()),
  d_leadDim(ref.LeadingDimension()),
  d_val_p(),
  d_work_p(),
  d_lwork(0)
  {
    
    size_t mylen = d_leadDim * d_numCol;
    d_val_p = new Scalar[mylen];
    d_hasAllocatedVal = true;
    for (size_t ii = 0; ii < mylen; ++ii)
      d_val_p[ii] = ref.d_val_p[ii];
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class T>
  DenseMatrix<Scalar>::DenseMatrix
  (
   const DenseMatrix<T> &ref
   )
  :
  d_numRow(ref.NumRows()),
  d_numCol(ref.NumVectors()),
  d_leadDim(ref.LeadingDimension()),
  d_val_p(),
  d_work_p(),
  d_lwork(0)
  {
    
    size_t mylen = d_leadDim * d_numCol;
    d_val_p = new Scalar[mylen];
    d_hasAllocatedVal = true;
    const T* rval_p = ref.Values();
    for (size_t ii = 0; ii < mylen; ++ii)
      d_val_p[ii] = (Scalar) rval_p[ii];
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>::~DenseMatrix
  (
  )
  {
    if ((d_hasAllocatedVal == true) && (d_val_p))
      delete[] d_val_p;
    d_val_p = (Scalar *) 0;
    
    if ((d_lwork > 0) && (d_work_p))
      delete[] d_work_p;
    d_work_p = (Scalar *) 0;
    d_lwork = 0;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class T>
  bool DenseMatrix<Scalar>::HasMatchingDimensions
  (
   const DenseMatrix<T> &xx
   ) const
  { 
    
    bool test = true;
    if ((d_numRow != xx.NumRows()) || (d_numCol != xx.NumVectors()))
      test = false;
    return test;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class Index>
  const DenseMatrix<Scalar> DenseMatrix<Scalar>::operator()
  (
   const std::vector<Index> &RowIndex,
   const std::vector<Index> &ColIndex
   )
  const
  {
    
    const size_t rsize = RowIndex.size();
    const size_t csize = ColIndex.size();
    
    //--- Check the row indices
    bool RowOutofBounds = false;
    for (size_t ir = 0; ir < rsize; ++ir)
    {
      if ((RowIndex[ir] >= d_numRow) || (RowIndex[ir] < 0))
      {
        RowOutofBounds = true;
        break;
      }
    }
    assert(RowOutofBounds == false);
    
    //--- Check the column indices
    bool ColOutofBounds = false;
    for (size_t ir = 0; ir < csize; ++ir)
    {
      if ((ColIndex[ir] >= d_numCol) || (ColIndex[ir] < 0))
      {
        ColOutofBounds = true;
        break;
      }
    }
    assert(ColOutofBounds == false);
    
    DenseMatrix<Scalar> res(rsize, csize);
    for (size_t ir = 0; ir < rsize; ++ir)
    {
      
      size_t ii = (size_t) RowIndex[ir];
      for (size_t jv = 0; jv < csize; ++jv)
      {
        size_t jj = (size_t) ColIndex[jv];
        res.d_val_p[ir + jv*res.d_numRow] = d_val_p[ii + jj*d_leadDim];
      }
      
    } // for (size_t ir = 0; ir < RowIndex.size(); ++ir)
    
    return res;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator=
  (
   const DenseMatrix<Scalar> &bv
   )
  {
    
    const Scalar *bval_p = bv.Values();
    size_t blen = bv.LeadingDimension() * bv.NumColumns();
    if ((d_numRow != bv.NumRows()) || 
        (d_numCol != bv.NumColumns()))
    {
      // Note that the leading dimensions can differ.
      // This choice allows to write in a sub-block of memory.
      if (d_hasAllocatedVal == true)
        delete[] d_val_p;
      d_val_p = new Scalar[blen];
      d_hasAllocatedVal = true;
      d_numRow = bv.NumRows();
      d_numCol = bv.NumColumns();
      d_leadDim = bv.LeadingDimension();
    }
    
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        d_val_p[ii + jj * d_leadDim] = bv(ii, jj);
      }
    }
    
    return *this;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class T>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator=
  (
   const DenseMatrix<T> &bv
   )
  {
    
    const T *bval_p = bv.Values();
    size_t blen = bv.LeadingDimension() * bv.NumColumns();
    if ((d_numRow != bv.NumRows()) || 
        (d_numCol != bv.NumColumns()))
    {
      // Note that the leading dimensions can differ.
      // This choice allows to write in a sub-block of memory.
      if (d_hasAllocatedVal == true)
        delete[] d_val_p;
      d_val_p = new Scalar[blen];
      d_hasAllocatedVal = true;
      d_numRow = bv.NumRows();
      d_numCol = bv.NumColumns();
      d_leadDim = bv.LeadingDimension();
    }
    
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        d_val_p[ii + jj * d_leadDim] = ((Scalar) bv(ii,jj));
      }
    }
    
    return *this;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator=
  (
   const Scalar &b
   )
  {
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        d_val_p[pos] = b;
      }
    }
    return *this;
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator*=
  (
   const Scalar &b
   )
  {
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        d_val_p[pos] *= b;
      }
    }
    return *this;
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator/=
  (
   const Scalar &b
   )
  {
    Scalar invb = ((Scalar) 1) / b;
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        d_val_p[pos] *= invb;
      }
    }
    return *this;
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator+=
  (
   const Scalar &b
   )
  {
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        d_val_p[pos] += b;
      }
    }
    return *this;
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class T>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator+=
  (
   const DenseMatrix<T> &bv
   )
  {
    
    bool MatchDim = this->HasMatchingDimensions(bv);
    if (MatchDim == false)
    {
      std::cerr << "\n !!! DenseMatrix += >> ";
      std::cerr << "the dimensions of vectors do not match";
      std::cerr << "!!! \n\n";
      assert(MatchDim == true);
    }
    
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        d_val_p[ii + jj * d_leadDim] += ((Scalar) bv(ii,jj));
      }
    }
    
    return *this;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator-=
  (
   const Scalar &b
   )
  {
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        d_val_p[pos] -= b;
      }
    }
    return *this;
  }
  
  
  //----------------------------------------------------------
  
  
  template <class Scalar>
  template <class T>
  DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator-=
  (
   const DenseMatrix<T> &bv
   )
  {
    
    bool MatchDim = this->HasMatchingDimensions(bv);
    if (MatchDim == false)
    {
      std::cerr << "\n !!! DenseMatrix -= >> ";
      std::cerr << "the dimensions of vectors do not match";
      std::cerr << "!!! \n\n";
      assert(MatchDim == true);
    }
    
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        d_val_p[ii + jj * d_leadDim] -= ((Scalar) bv(ii,jj));
      }
    }
    
    return *this;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  int DenseMatrix<Scalar>::Add
  (
   int gRow,
   int gCol,
   Scalar val
   )
  {
    
    if ((gRow < 0) || (gRow >= d_numRow))
      return -1; 
    
    if ((gCol < 0) || (gRow >= d_numCol))
      return -1; 
    
    d_val_p[gRow + gCol * d_leadDim] += val;
    
    return 0;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  int DenseMatrix<Scalar>::Add
  (
   const std::vector<int> &pos,
   const std::vector<Scalar> &squareMat
   )
  {
    
    bool isGood = true;
    
    int iLen = pos.size();
    if (squareMat.size() != iLen * iLen)
      return -1;
    
    for (size_t ii = 0; ii < iLen; ++ii)
    {
      if ((pos[ii] < 0) || (pos[ii] >= d_numRow) || (pos[ii] >= d_numCol))
      {
        isGood = false;
        break;
      }
    }
    
    if (isGood == false)
      return -1;
    
    for (int ii = 0; ii < iLen; ++ii)
    {
      int iRow = pos[ii];
      for (int jj = 0; jj < iLen; ++jj)
      {
        int jCol = pos[jj];
        d_val_p[iRow + jCol * d_leadDim] += squareMat[ii + jj * iLen];
      }
    }
    
    return 0;
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  int DenseMatrix<Scalar>::Apply
  (
   const BlockVector<Scalar> &X,
   BlockVector<Scalar> &Y
   )
  {
   
    if ((Y.NumRows() != d_numRow) || (X.NumRows() != d_numCol)
        || (Y.NumVectors() < X.NumVectors()))
      return -1;
    
    Scalar *xval = X.Values();
    Scalar *yval = Y.Values();
    size_t xcol = X.NumVectors();
    size_t xrow = X.NumRows();
    
    for (size_t ii = 0; ii < d_numRow; ++ii)
    {
      for (size_t jj = 0; jj < xcol; ++jj)
      {
        Scalar tmp = 0;
        for (size_t kk = 0; kk < d_numCol; ++kk)
        {
          tmp += d_val_p[ii + kk * d_leadDim] * xval[kk + jj * xrow];
        }
        yval[ii + jj * d_numRow] = tmp;
      }
    }
    
    return 0;
    
  }
  
  
  // This line is essential to declare the specialization.
  // Without it, the code will compile but it will not use the specialized version.
  template<> int DenseMatrix<double>::Apply
  ( const BlockVector<double> &X, BlockVector<double> &Y);
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  template<class T>
  void DenseMatrix<Scalar>::AXPBY
  (
   const Scalar &a,
   const DenseMatrix<T> &X,
   const Scalar &b
   )
  {
    
    bool MatchDim = this->HasMatchingDimensions(X);
    if (MatchDim == false)
    {
      std::cerr << "\n !!! DenseMatrix AXPBY >> ";
      std::cerr << "the dimensions of vectors do not match";
      std::cerr << "!!! \n\n";
      assert(MatchDim == true);
    }
    
    if (b == (Scalar) 0)
    {
      if (a == (Scalar) 0)
      {
        for (size_t jj = 0; jj < d_numCol; ++jj)
        {
          for (size_t ii = 0; ii < d_numRow; ++ii)
          {
            d_val_p[ii + jj * d_leadDim] = a;
          }
        }
      }
      else
      {
        for (size_t jj = 0; jj < d_numCol; ++jj)
        {
          for (size_t ii = 0; ii < d_numRow; ++ii)
          {
            d_val_p[ii + jj * d_leadDim] = a * ((Scalar) X(ii,jj));
          }
        }
      }
    }
    else if (a == (Scalar) 0)
    {
      for (size_t jj = 0; jj < d_numCol; ++jj)
      {
        for (size_t ii = 0; ii < d_numRow; ++ii)
        {
          d_val_p[ii + jj * d_leadDim] *= b;
        }
      }
    }
    else
    {
      for (size_t jj = 0; jj < d_numCol; ++jj)
      {
        for (size_t ii = 0; ii < d_numRow; ++ii)
        {
          Scalar &myval = d_val_p[ii + jj * d_leadDim];
          myval = a * ((Scalar) X(ii,jj)) + b * myval;
        }
      }
    }
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  template<class T, class Index>
  void DenseMatrix<Scalar>::local_AXPBY
  (
   const std::vector<Index> &RowIndex,
   const Scalar &a,
   const DenseMatrix<T> &X,
   const Scalar &b
   )
  {
    
    const size_t xrow = X.NumRows();
    const size_t xcol = X.NumCols();
    const size_t rsize = RowIndex.size();
    
    if ((rsize != xrow) || (d_numCol != xcol))
    {
      bool MatchDim = false;
      assert(MatchDim == true); 
    }
    
    //--- Check the row indices
    bool RowOutofBounds = false;
    for (size_t ir = 0; ir < rsize; ++ir)
    {
      if ((RowIndex[ir] >= d_numRow) || (RowIndex[ir] < 0))
      {
        RowOutofBounds = true;
        break;
      } 
    }
    assert(RowOutofBounds == false);
    
    if (b == 0)
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < d_numCol; ++jc)
        {
          size_t pos = ii + jc*d_leadDim;
          d_val_p[pos] += a * ((Scalar) X(ir,jc));
        }
      }
    }
    else if (a == 0)
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < d_numCol; ++jc)
        {
          size_t pos = ii + jc*d_leadDim;
          d_val_p[pos] *= b;
        }
      }
    }
    else
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < d_numCol; ++jc)
        {
          size_t pos = ii + jc * d_leadDim;
          d_val_p[pos] *= b;
          d_val_p[pos] += a * ((Scalar) X(ir,jc));
        }
      }
    }
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  template<class T, class Index>
  void DenseMatrix<Scalar>::local_AXPBY
  (
   const std::vector<Index> &RowIndex,
   const std::vector<Index> &ColIndex,
   const Scalar &a,
   const DenseMatrix<T> &X,
   const Scalar &b
   )
  {
    
    const size_t xrow = X.NumRows();
    const size_t xcol = X.NumColumns();
    const size_t rsize = RowIndex.size();
    const size_t csize = ColIndex.size();
    
    if ((rsize != xrow) || (csize != xcol))
    {
      bool MatchDim = false;
      assert(MatchDim == true); 
    }
    
    //--- Check the row indices
    bool RowOutofBounds = false;
    for (size_t ir = 0; ir < rsize; ++ir)
    {
      if ((RowIndex[ir] >= d_numRow) || (RowIndex[ir] < 0))
      {
        RowOutofBounds = true;
        break;
      } 
    }
    assert(RowOutofBounds == false);
    
    //--- Check the column indices
    bool ColOutofBounds = false;
    for (size_t ir = 0; ir < csize; ++ir)
    {
      if ((ColIndex[ir] >= d_numCol) || (ColIndex[ir] < 0))
      {
        ColOutofBounds = true;
        break;
      } 
    }
    assert(ColOutofBounds == false);
    
    if (b == 0)
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < csize; ++jc)
        {
          size_t pos = ii + ((size_t) ColIndex[jc]) * d_leadDim;
          d_val_p[pos] += a * ((Scalar) X(ir,jc));
        }
      }
    }
    else if (a == 0)
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < csize; ++jc)
        {
          size_t pos = ii + ((size_t) ColIndex[jc]) * d_leadDim;
          d_val_p[pos] *= b;
        }
      }
    }
    else
    {
      for (size_t ir = 0; ir < rsize; ++ir)
      {
        size_t ii = (size_t) RowIndex[ir];
        for (size_t jc = 0; jc < csize; ++jc)
        {
          size_t pos = ii + ((size_t) ColIndex[jc]) * d_leadDim;
          d_val_p[pos] *= b;
          d_val_p[pos] += a * ((Scalar) X(ir,jc));
        }
      }
    }
    
  }
  
  
  //----------------------------------------------------------
  
  
  template<class Scalar>
  void DenseMatrix<Scalar>::Random
  (
   int seed
   ) 
  {
    
    if (seed == 0) 
    {
      if (MathCXX::d_hasInitializedSRand == false) 
      {
        srand((unsigned int) time(0));
        MathCXX::d_hasInitializedSRand = true;
      } 
    }
    else if (seed < 0) 
    {
      srand((unsigned int) time(0));
    }
    else 
    {
      srand((unsigned int) seed);
    }
    
    size_t pos;
    for (size_t jj = 0; jj < d_numCol; ++jj)
    {
      for (size_t ii = 0; ii < d_numRow; ++ii)
      {
        pos = ii + jj * d_leadDim;
        SetToRandom(d_val_p[pos]);
      }
    }
    
  }
  
  
  //----------------------------------------------------------
  
}


#endif
