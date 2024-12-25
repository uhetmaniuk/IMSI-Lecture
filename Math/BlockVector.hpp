#ifndef MATH_CXX_BLOCKVECTOR_HPP
#define MATH_CXX_BLOCKVECTOR_HPP 


#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <sys/time.h>
#include <vector>


//class Communicator;


namespace MathCXX {
  
  
  template<class Scalar> class DenseMatrix;
  
  
  //! A templated class for block vectors.
  
  /// The BlockVector class enables the use of block vectors.
  template <class Scalar>
  class BlockVector {
    
  private:
    
    
  protected:
    
    DenseMatrix<Scalar> *d_val_p;
    
    
    //--- Data for the distribution
//    bool d_isDistributed;
//    const Communicator *d_Comm_p;
//#ifdef FEC_USE_TRILINOS
//    Epetra_BlockMap *d_RowMap_p;
//#endif
    
    
  public:
    
    
    ///////////
    // CREATOR
    ///////////
    
    
    //! Constructor.
    //! \param[in] nRow Number of rows for each vector.
    //! \param[in] nVec Number of vectors in the block.
    /// \note The values of the multivector are not initialized.
    BlockVector<Scalar>
    (
     size_t nRow, 
     size_t nVec
     );
    
    
    //! Destructor.
    virtual ~BlockVector
    (
    );
    
    
    ///////////////
    // MANIPULATOR 
    ///////////////
    
    
    //! Set all values in a multi-vector with constant value.
    /// \param[in] arg Value to use.
    /// \return Updated multi-vector.
    BlockVector<Scalar>& operator=
    (
     const Scalar &arg
     );
    
    
    //! Sets block vector values to random numbers.
    /// \param[in] seed Integer setting the state of the random generator.
    /// \note A seed negative or zero will reset the state to the time.
    /// \note The real entries will be between -1 and 1.
    void Random
    (
      int seed = 0
    );
    
    
    ////////////
    // ACCESSOR 
    ////////////
    
    
    //! Returns the number of rows in the block.
    size_t NumRows
    () const;
    
    
    //! Returns the number of vectors in the block.
    size_t NumVectors
    () const;
    
    
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
    
    
    //! Returns a reference to the pointer of values.
    /// \note This function gives a read-only access to the pointer.
    Scalar* Values
    () const;
    
    
    //////////
    // FRIEND
    //////////

    
  };
  
  
}


//==============================================================================
//==============================================================================



#include "DenseMatrix.hpp"


namespace MathCXX 
{
  
  
  /////////////////////////////////////
  //
  // Definition of templated member functions
  //
  ////////////////////////////////////
  
  
  
  //------------------------------------------------------------------
  
  template<class Scalar>
  BlockVector<Scalar>::BlockVector
  (
   size_t nRow, 
   size_t nVec
   )
  : 
  d_val_p(new DenseMatrix<Scalar>(nRow, nVec))
  //  , d_isDistributed(false)
  //  , d_Comm_p(0)
  //#ifdef FEC_USE_TRILINOS
  //  , d_RowMap_p(0)
  //#endif
  {
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  BlockVector<Scalar>::~BlockVector
  (
  )
  {
    
    delete d_val_p;
    d_val_p = (DenseMatrix<Scalar> *) 0;
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  BlockVector<Scalar>& BlockVector<Scalar>::operator=
  (
   const Scalar &b
   )
  {
    
    *d_val_p = b;
    
    return *this;
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  void BlockVector<Scalar>::Random
  (
   int seed
   )
  {
    
    d_val_p->Random(seed);
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  size_t BlockVector<Scalar>::NumRows
  (
  )
  const
  {
    
    return d_val_p->NumRows();
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  size_t BlockVector<Scalar>::NumVectors
  (
  )
  const
  {
    
    return d_val_p->NumColumns();
    
  } 
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  Scalar& BlockVector<Scalar>::operator() 
  (
   size_t RowIndex, 
   size_t ColIndex
   )
  {
    
    return (*d_val_p)(RowIndex, ColIndex);
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  Scalar BlockVector<Scalar>::operator() 
  (
   size_t RowIndex, 
   size_t ColIndex
   )
  const
  {
    
    return (*d_val_p)(RowIndex, ColIndex);
    
  }
  
  //------------------------------------------------------------------
  
  template <class Scalar>
  Scalar* BlockVector<Scalar>::Values
  (
  )
  const
  {
    
    return d_val_p->Values();
    
  }
  
  //------------------------------------------------------------------
  
  
}



#endif
