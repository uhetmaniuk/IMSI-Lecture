#pragma once

#include <complex>
#include <iostream>
#include <Kokkos_SIMD.hpp>

//----------------------------------------------------
//
// Class Definition
//
//----------------------------------------------------

template <class Scalar>
class SparseMatrix {
  
protected: 
  
  int d_numRow;
  int d_numCol;
  int d_NNZ;
  
  int *d_ROWPTR;
  int *d_COLIDX;
  Scalar *d_ANZ;
  
  bool d_hasAllocatedData;
  
public:
  
  //////////////////////////////
  // Creator
  //////////////////////////////
  SparseMatrix
  (
   int nrow,
   int ncol,
   int NNZ, 
   int ROWPTR[], 
   int COLIDX[],
   Scalar ANZ[],
   bool CopyData = false
   );
  
  SparseMatrix
  (
   int nrow,
   int ncol,
   int NNZ
   );
  
  ~SparseMatrix();
  
  //////////////////////////////
  // Manipulator
  //////////////////////////////
  
  void cleanup();
  
  int Apply
  (
   int NRHS,
   Scalar X[],
   Scalar AX[]
   )  const;
  
  //! 
  //! ...
  //! This routine extracts a submatrix with same indices on rows and columns.
  //! The resulting submatrix will be symmetric and sparse.
  //! The arrays for indices and data are copied.
  //!
  SparseMatrix<Scalar>* Extract
  (
   int row_nIndex,
   int row_INDEX[],
   int col_nIndex,
   int col_INDEX[]
   ) const;
  
  
  int ReplaceValue
  (
      int row,
      int col,
      Scalar val
  ) const; 

  
  //////////////////////////////
  // Accessor
  //////////////////////////////
  
  int* GetRowPointers() const { return d_ROWPTR; };
  
  int* GetColIndices() const { return d_COLIDX; };
  
  Scalar* GetData() const { return d_ANZ; };
  
  double GetMaxNorm() const;
  
  int GetNumCol() const { return d_numCol; };
  
  int GetNumNonZeros() const { return d_NNZ; };

  int GetNumRow() const { return d_numRow; };
  
  
  //////////////////////////////
  // Friend functions
  //////////////////////////////
  
  
  template <class MyType> friend MyType* Full(SparseMatrix<MyType> &K);
  
};


//-------------------------------------------------------------------
//
// Definition of the functions for the templates class
//
//-------------------------------------------------------------------


#include <cmath>


template <class Scalar>
SparseMatrix<Scalar>::SparseMatrix
(
 int nrow,
 int ncol,
 int NNZ, 
 int ROWPTR[], 
 int COLIDX[],
 Scalar ANZ[],
 bool CopyData
)
{
  
  //
  // Input:
  //   nrow = number of rows
  //   ncol = number of columns
  //   NNZ = number of nonzeros in full matrix
  //   COLIDX[ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero column numbers for row i
  //   ANZ[   ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero values in row i
  //   Note: C-style numbering of inputs is assumed 
  //   (i.e. row and column numbers are all between 0 and N_-1)
  //
  
  d_numRow = nrow;
  d_numCol = ncol;
  d_NNZ = NNZ;
  
  if (CopyData == false)
  {
    d_ROWPTR = ROWPTR;
    d_COLIDX = COLIDX;
    d_ANZ = ANZ;
    d_hasAllocatedData = false;
  }
  else
  {
    d_ROWPTR = new int[d_numRow+1];
    for (int jj = 0; jj <= d_numRow; ++jj)
    {
      d_ROWPTR[jj] = ROWPTR[jj];
    }
    d_COLIDX = new int[d_NNZ];
    d_ANZ = new Scalar[d_NNZ];
    for (int ii = 0; ii < d_NNZ; ++ii)
    {
      d_COLIDX[ii] = COLIDX[ii];
      d_ANZ[ii] = ANZ[ii];
    }
    d_hasAllocatedData = true;
  }
  
}


template <class Scalar>
SparseMatrix<Scalar>::SparseMatrix
(
 int nrow,
 int ncol,
 int NNZ
)
{
  
  //
  // Input:
  //   nrow = number of rows
  //   ncol = number of columns
  //   NNZ = number of nonzeros in full matrix
  //   COLIDX[ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero column numbers for row i
  //   ANZ[   ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero values in row i
  //   Note: C-style numbering of inputs is assumed 
  //   (i.e. row and column numbers are all between 0 and N_-1)
  //
  
  d_numRow = nrow;
  d_numCol = ncol;
  d_NNZ = NNZ;
  
  d_ROWPTR = new int[d_numRow+1];
  for (int ii = 0; ii <= d_numRow; ++ii)
    d_ROWPTR[ii] = 0;
  
  d_COLIDX = new int[d_NNZ];
  d_ANZ = new Scalar[d_NNZ];
  for (int ii = 0; ii < d_NNZ; ++ii)
  {
    d_COLIDX[ii] = -1;
    d_ANZ[ii] = (Scalar) 0;
  }
  
  d_hasAllocatedData = true;
  
}


template <class Scalar>
SparseMatrix<Scalar>::~SparseMatrix
(
)
{
  
  if (d_hasAllocatedData == true)
  {
    if (d_ROWPTR)
      delete[] d_ROWPTR;
    d_ROWPTR = 0;
    if (d_COLIDX)
      delete[] d_COLIDX;
    d_COLIDX = 0;
    if (d_ANZ)
      delete[] d_ANZ;
    d_ANZ = 0;
  } //  if (d_hasAllocatedData == true)
  
}


template <class Scalar>
void SparseMatrix<Scalar>::cleanup
(
)
{
  
  if (d_hasAllocatedData == true)
  {
    if (d_ROWPTR)
      delete[] d_ROWPTR;
    d_ROWPTR = 0;
    if (d_COLIDX)
      delete[] d_COLIDX;
    d_COLIDX = 0;
    if (d_ANZ)
      delete[] d_ANZ;
    d_ANZ = 0;
  } //  if (d_hasAllocatedData == true)
  
}


template <class Scalar>
int SparseMatrix<Scalar>::Apply
(
 int NRHS,
 Scalar X[],
 Scalar AX[]
)  const
{
  using simd_t        = Kokkos::Experimental::native_simd<Scalar>;
  constexpr int width = simd_t::size();

  for (int ii = 0; ii < d_numRow; ++ii) {
    int kk = 0;
    for (; kk <= NRHS - width; kk += width) {
      simd_t sum(0);
      for (int jj = d_ROWPTR[ii]; jj < d_ROWPTR[ii + 1]; ++jj) {
        Scalar val   = d_ANZ[jj];
        int    coljj = d_COLIDX[jj];
        simd_t x_v;
        for (int l = 0; l < width; ++l) { x_v[l] = X[coljj + (kk + l) * d_numRow]; }
        sum += val * x_v;
      }
      for (int l = 0; l < width; ++l) { AX[ii + (kk + l) * d_numRow] = sum[l]; }
    }
    for (; kk < NRHS; ++kk) {
      Scalar sum(0);
      for (int jj = d_ROWPTR[ii]; jj < d_ROWPTR[ii + 1]; ++jj) {
        Scalar val   = d_ANZ[jj];
        int    coljj = d_COLIDX[jj];
        sum += val * X[coljj + kk * d_numRow];
      }
      AX[ii + kk * d_numRow] = sum;
    }
  }  // for (int ii = 0; ii < d_numRow; ++ii)
  return 0;
}


template <class Scalar>
SparseMatrix<Scalar>* SparseMatrix<Scalar>::Extract
(
 int row_nIndex,
 int row_INDEX[],
 int col_nIndex,
 int col_INDEX[]
 ) const
{
  
  std::vector<int> row_globalToKeep(d_numRow, -1);
  for (int ii = 0; ii < row_nIndex; ++ii)
  {
    row_globalToKeep[row_INDEX[ii]] = ii;
  }
  
  std::vector<int> col_globalToKeep(d_numCol, -1);
  for (int ii = 0; ii < col_nIndex; ++ii)
  {
    col_globalToKeep[col_INDEX[ii]] = ii;
  }
  
  int my_nnz = 0;
  for (int ii = 0; ii < d_numRow; ++ii)
  {
    if (row_globalToKeep[ii] == -1)
      continue;
    for (int jj = d_ROWPTR[ii]; jj < d_ROWPTR[ii+1]; ++jj)
    {
      if (col_globalToKeep[d_COLIDX[jj]] > -1)
        my_nnz += 1;
    }
  }
  
  SparseMatrix<Scalar> *MyK = new SparseMatrix<Scalar>(row_nIndex, col_nIndex, my_nnz);
  
  int *my_rowbeg = MyK->GetRowPointers();
  int *my_colidx = MyK->GetColIndices();
  Scalar *my_data = MyK->GetData();
  
  int count = 0;
  my_nnz = 0;
  for (int ii = 0; ii < d_numRow; ++ii)
  {
    if (row_globalToKeep[ii] == -1)
      continue;
    my_rowbeg[count+1] = my_rowbeg[count];
    for (int jj = d_ROWPTR[ii]; jj < d_ROWPTR[ii+1]; ++jj)
    {
      if (col_globalToKeep[d_COLIDX[jj]] > -1)
      {
        my_rowbeg[count+1] += 1;
        my_colidx[my_nnz] = col_globalToKeep[d_COLIDX[jj]];
        my_data[my_nnz] = d_ANZ[jj];
        my_nnz += 1;
      }
    }
    count += 1;
  }
  
  return MyK;
  
}


template <class Scalar>
int SparseMatrix<Scalar>::ReplaceValue
(
  int row,
  int col,
  Scalar val
)
const
{
  
  int info = 0;
  
  if ((row < 0) || (row >= d_numRow))
    return -1;
  
  if ((col < 0) || (col >= d_numCol))
    return -1;

  info = -2;
  for (int jj = d_ROWPTR[row]; jj < d_ROWPTR[row+1]; ++jj)
  {
    if ((d_COLIDX[jj] == col) || (d_COLIDX[jj] == -1))
    {
      d_COLIDX[jj] = col;
      d_ANZ[jj] = val;
      info = (d_COLIDX[jj] == -1) ? 1 : 0;
      break;
    }
  }
  
  return info;
  
}



template <class Scalar>
double SparseMatrix<Scalar>::GetMaxNorm
(
)
const
{
  int i, j;
  double anorm = 0;
  for (i = 0; i < d_numRow; i++) {
    double t = 0;
    for (j = d_ROWPTR[i]; j < d_ROWPTR[i+1]; j++) {
      t += std::abs(d_ANZ[j]);
    }
    if (t > anorm)
      anorm=t;
  }
  return anorm;
}


  
template <class MyType>
MyType* Full
(
 SparseMatrix<MyType> &K
)
{
  
  int numRows = K.GetNumRow();
  int numCols = K.GetNumCol();
  
  int numEntries = numRows * numCols;
  MyType *KK = new MyType[numEntries];
  for (int ii = 0; ii < numEntries; ++ii)
    KK[ii] = (MyType) 0;

  int *rowbeg = K.GetRowPointers();
  int *colidx = K.GetColIndices();
  MyType *data = K.GetData();

  for (int ii = 0; ii < numRows; ++ii)
  {
    for (int jj = rowbeg[ii]; jj < rowbeg[ii + 1]; ++jj)
    {
      int mycol = colidx[jj];
      KK[ii + mycol * numRows] = data[jj];
    }
  }
  
  return KK;
  
}


