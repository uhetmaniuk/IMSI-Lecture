#pragma once

#include <complex>
#include <iostream>
#include <set>
#include <vector>

#include "SparseMatrix.hpp"
#include "my_sparse_solver.h"

extern "C" {

void
metis_nodend(int* N, int XADJ[], int ADJD[], int* numflag, int OPTIONS[], int PERM[], int INVP[]);
}

//----------------------------------------------------
//
// Class Definition
//
//----------------------------------------------------

template <class Scalar>
class SymmetricSparse : public SparseMatrix<Scalar>
{
 protected:
  //---------------
  template <class MyScalar>
  class SparseCholeskyData
  {
   public:
    int       d_DEFBLK, d_NSUPER, d_NDEF, d_LBDEF;
    int*      d_XSUPER;
    int*      d_XLINDX;
    int*      d_LINDX;
    int*      d_XLNZ;
    int*      d_PERM;
    int*      d_INVP;
    int*      d_IPROW;
    int*      d_IPCOL;
    int*      d_DEF;
    MyScalar* d_LNZ;
    MyScalar* d_NS;

    SparseCholeskyData()
    {
      d_DEFBLK = 0;
      d_NSUPER = 0;
      d_NDEF   = 0;
      d_LBDEF  = 0;
      d_XSUPER = 0;
      d_XLINDX = 0;
      d_LINDX  = 0;
      d_XLNZ   = 0;
      d_PERM   = 0;
      d_INVP   = 0;
      d_IPROW  = 0;
      d_IPCOL  = 0;
      d_DEF    = 0;
      //---
      d_LNZ = 0;
      d_NS  = 0;
    };
  };
  //---------------

  SparseCholeskyData<Scalar>* d_Factor;

  /////////////////////////////////
  // Protected Manipulator
  /////////////////////////////////

  void
  inpnv(
      int&   n,
      int    colptr[],
      int    rowidx[],
      Scalar values[],
      int    perm[],
      int    invp[],
      int&   nsuper,
      int    xsuper[],
      int    xlindx[],
      int    lindx[],
      int    xlnz[],
      Scalar lnz[],
      int    offset[]);

  void
  numericalFactorization(int*& SNODE, int*& IWORK, int& IWSIZE, int& IFLAG);

  void
  symbolicFactorization(int*& SNODE, int*& IWORK, int& IWSIZE, int& IFLAG);

 public:
  //////////////////////////////
  // Creator
  //////////////////////////////
  SymmetricSparse(int N_, int NNZ, int ROWPTR[], int COLIDX[], Scalar ANZ[], bool CopyData = false);

  SymmetricSparse(int N_, int NNZ);

  ~SymmetricSparse();

  //////////////////////////////
  // Manipulator
  //////////////////////////////
  void
  cleanup();

  int
  factor();

  int
  Solve(int NRHS, Scalar RHS[], Scalar SOL[], Scalar TEMP[]);

  int
  Solve(int NRHS, Scalar RHS[], Scalar SOL[]);

  //!
  //! ...
  //! This routine extracts a submatrix with same indices on rows and columns.
  //! The resulting submatrix will be symmetric and sparse.
  //! The arrays for indices and data are copied.
  //!
  SymmetricSparse<Scalar>*
  SymmetricExtract(int nIndex, int INDEX[]) const;

  //!
  //! ...
  //!
  SymmetricSparse<Scalar>*
  ExtractSchurComplement(std::vector<int>& keepRows, std::vector<int>& elimRows) const;

  //!
  //! ...
  //!
  SymmetricSparse<Scalar>*
  ExtractSchurComplement(int numKeep, int keepRows[], int numElim, int elimRows[]) const;

  //////////////////////////////
  // Accessor
  //////////////////////////////
};

//-------------------------------------------------------------------
//
// Definition of the functions for the templates class
//
//-------------------------------------------------------------------

#include <cmath>

template <class Scalar>
SymmetricSparse<Scalar>::SymmetricSparse(int N_, int NNZ, int ROWPTR[], int COLIDX[], Scalar ANZ[], bool CopyData)
    : SparseMatrix<Scalar>(N_, N_, NNZ, ROWPTR, COLIDX, ANZ, CopyData)
{
  /*
   //
   // Input:
   //   N_ = number of equations
   //   NNZ = number of nonzeros in full matrix
   //   COLIDX[ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero column numbers for row i
   //   ANZ[   ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero values in row i
   //   Note: C-style numbering of inputs is assumed
   //   (i.e. row and column numbers are all between 0 and N_-1)
   //   Note: For symmetric matrices, (ROWPTR, COLIDX) and (COLPTR, ROWIDX)
   //   are the same.
   //

   d_N = N_;
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
   d_ROWPTR = new int[d_N+1];
   for (int jj = 0; jj <= d_N; ++jj)
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
   */

  d_Factor = 0;
}

template <class Scalar>
SymmetricSparse<Scalar>::SymmetricSparse(int N_, int NNZ) : SparseMatrix<Scalar>(N_, N_, NNZ)
{
  /*
   //
   // Input:
   //   N_ = number of equations
   //   NNZ = number of nonzeros in full matrix
   //   COLIDX[ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero column numbers for row i
   //   ANZ[   ROWPTR[i]:ROWPTR[i+1]-1]  = nonzero values in row i
   //   Note: C-style numbering of inputs is assumed
   //   (i.e. row and column numbers are all between 0 and N_-1)
   //   Note: For symmetric matrices, (ROWPTR, COLIDX) and (COLPTR, ROWIDX)
   //   are the same.
   //

   d_N = N_;
   d_NNZ = NNZ;

   d_ROWPTR = new int[d_N+1];
   for (int ii = 0; ii <= d_N; ++ii)
   d_ROWPTR[ii] = 0;

   d_COLIDX = new int[d_NNZ];
   d_ANZ = new Scalar[d_NNZ];
   for (int ii = 0; ii < d_NNZ; ++ii)
   {
   d_COLIDX[ii] = 0;
   d_ANZ[ii] = (Scalar) 0;
   }

   d_hasAllocatedData = true;
   */

  d_Factor = 0;
}

template <class Scalar>
SymmetricSparse<Scalar>::~SymmetricSparse()
{
  if (d_Factor) {
    if (d_Factor->d_XSUPER) delete[] d_Factor->d_XSUPER;
    d_Factor->d_XSUPER = 0;
    if (d_Factor->d_XLINDX) delete[] d_Factor->d_XLINDX;
    d_Factor->d_XLINDX = 0;
    if (d_Factor->d_LINDX) delete[] d_Factor->d_LINDX;
    d_Factor->d_LINDX = 0;
    if (d_Factor->d_XLNZ) delete[] d_Factor->d_XLNZ;
    d_Factor->d_XLNZ = 0;
    if (d_Factor->d_PERM) delete[] d_Factor->d_PERM;
    d_Factor->d_PERM = 0;
    if (d_Factor->d_INVP) delete[] d_Factor->d_INVP;
    d_Factor->d_INVP = 0;
    if (d_Factor->d_IPROW) delete[] d_Factor->d_IPROW;
    d_Factor->d_IPROW = 0;
    if (d_Factor->d_IPCOL) delete[] d_Factor->d_IPCOL;
    d_Factor->d_IPCOL = 0;
    if (d_Factor->d_DEF) delete[] d_Factor->d_DEF;
    d_Factor->d_DEF = 0;
    if (d_Factor->d_LNZ) delete[] d_Factor->d_LNZ;
    d_Factor->d_LNZ = 0;
    if (d_Factor->d_NS) delete[] d_Factor->d_NS;
    d_Factor->d_NS = 0;
  }
  delete d_Factor;
  d_Factor = 0;
}

template <class Scalar>
void
SymmetricSparse<Scalar>::cleanup()
{
  SparseMatrix<Scalar>::cleanup();

  if (d_Factor) {
    if (d_Factor->d_XSUPER) delete[] d_Factor->d_XSUPER;
    d_Factor->d_XSUPER = 0;
    if (d_Factor->d_XLINDX) delete[] d_Factor->d_XLINDX;
    d_Factor->d_XLINDX = 0;
    if (d_Factor->d_LINDX) delete[] d_Factor->d_LINDX;
    d_Factor->d_LINDX = 0;
    if (d_Factor->d_XLNZ) delete[] d_Factor->d_XLNZ;
    d_Factor->d_XLNZ = 0;
    if (d_Factor->d_PERM) delete[] d_Factor->d_PERM;
    d_Factor->d_PERM = 0;
    if (d_Factor->d_INVP) delete[] d_Factor->d_INVP;
    d_Factor->d_INVP = 0;
    if (d_Factor->d_IPROW) delete[] d_Factor->d_IPROW;
    d_Factor->d_IPROW = 0;
    if (d_Factor->d_IPCOL) delete[] d_Factor->d_IPCOL;
    d_Factor->d_IPCOL = 0;
    if (d_Factor->d_DEF) delete[] d_Factor->d_DEF;
    d_Factor->d_DEF = 0;
    if (d_Factor->d_LNZ) delete[] d_Factor->d_LNZ;
    d_Factor->d_LNZ = 0;
    if (d_Factor->d_NS) delete[] d_Factor->d_NS;
    d_Factor->d_NS = 0;
  }
  delete d_Factor;
  d_Factor = 0;
}

template <class Scalar>
int
SymmetricSparse<Scalar>::factor()
{
  int  IWSIZE, IFLAG;
  int* SNODE = new int[this->d_numRow];
  int* IWORK;

  //--- Create structure to store the Cholesky factorization
  d_Factor = new SparseCholeskyData<Scalar>();

  if (this->d_numRow == 1) {
    //--- UH (06/2009)
    // This is a fix when n = 1. There was a bug on the Mac.
    // The modification might not be consistent with the general case.
    //--------
    d_Factor->d_XLNZ    = new int[this->d_numRow + 1];
    d_Factor->d_XLNZ[0] = 0;
    d_Factor->d_XLNZ[1] = 1;
    int NLNZ            = d_Factor->d_XLNZ[this->d_numRow];
    d_Factor->d_LNZ     = new Scalar[NLNZ];
    d_Factor->d_LNZ[0]  = this->d_ANZ[0];
    if (this->d_ANZ[0] == (Scalar)0) {
      std::cout << "\n !!! SymmetricSparse<Scalar>::factor >> The matrix is singular !!! \n\n" << std::endl;
      return -1;
    }
    d_Factor->d_DEF = new int[this->d_numRow];
    for (int i = 0; i < this->d_numRow; i++) d_Factor->d_DEF[i] = 0;
    d_Factor->d_NDEF  = 0;
    d_Factor->d_LBDEF = 0;
    //-------
    return 0;
  }

  //
  //--- Perform the symbolic factorization
  //--- convert COLPTR, ROWIDX, XADJ, and ADJ to Fortran numbering
  //

  symbolicFactorization(SNODE, IWORK, IWSIZE, IFLAG);

  //--- Perform the numerical factorization
  numericalFactorization(SNODE, IWORK, IWSIZE, IFLAG);

  delete[] IWORK;
  delete[] SNODE;

  //
  // convert ROWPTR and COLIDX back to C numbering
  //

  for (int i = 0; i <= this->d_numRow; i++) this->d_ROWPTR[i]--;
  for (int i = 0; i < this->d_NNZ; i++) this->d_COLIDX[i]--;

  return 0;
}

template <class Scalar>
void
SymmetricSparse<Scalar>::numericalFactorization(int*& SNODE, int*& IWORK, int& IWSIZE, int& IFLAG)
{
  //
  //--- Start working with matrix values
  //
  // input numerical values into data structures
  //

  int TMPSIZ, RWSIZE;

  int NLNZ = d_Factor->d_XLNZ[this->d_numRow];
  //  std::cout << "number of nonzeros in LU factorization = " << NLNZ << std::endl;
  //  std::cout << "NLNZ = " << NLNZ << std::endl;
  d_Factor->d_LNZ = new Scalar[NLNZ];
  for (int i = 0; i < NLNZ; i++) { d_Factor->d_LNZ[i] = 0; }
  for (int i = 0; i < this->d_numRow; i++) { d_Factor->d_DEF[i] = 0; }
  d_Factor->d_NDEF  = 0;
  d_Factor->d_LBDEF = 0;

  //       ***************************************************
  //       Numerical input into data structure for sparse LDL'
  //       factorization.
  //       ***************************************************
  //
  //       --------------------------------------------------------
  //       INPNV ...   input numerical values into data structures.
  //
  //       Input:      N, COLPTR, ROWIDX, VALUES, PERM, INVP,
  //                   NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
  //       Output:     LNZ
  //       Work:       IWORK(N)
  //       --------------------------------------------------------

  //
  // Note (10/18/08 -- UH)
  // When calling invpnv, we exploit the symmetric property.
  // invpnv assumes that the matrix is stored column-wise.
  // For a symmetric matrix fully-stored (not the upper half or lower-half)
  // it does not matter which arrays are passed.
  //

  SymmetricSparse<Scalar>::inpnv(
      this->d_numRow,
      this->d_ROWPTR,
      this->d_COLIDX,
      this->d_ANZ,
      d_Factor->d_PERM,
      d_Factor->d_INVP,
      d_Factor->d_NSUPER,
      d_Factor->d_XSUPER,
      d_Factor->d_XLINDX,
      d_Factor->d_LINDX,
      d_Factor->d_XLNZ,
      d_Factor->d_LNZ,
      IWORK);

  //       ************************
  //       Numerical factorization.
  //       ************************
  //       ---------------------------------------------------
  //       BFINIT ...  initialization for block factorization.
  //
  //       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX
  //       Output:     TMPSIZ, RWSIZE
  //       ---------------------------------------------------

  bfinit_(d_Factor->d_NSUPER, d_Factor->d_XSUPER, SNODE, d_Factor->d_XLINDX, d_Factor->d_LINDX, TMPSIZ, RWSIZE);

  if (TMPSIZ < 1) { TMPSIZ = 1; }
  TMPSIZ           = 2 * TMPSIZ;
  Scalar* TMPVEC   = new Scalar[TMPSIZ];
  int     RWORKdim = this->d_numRow;
  if (RWSIZE > this->d_numRow) RWORKdim = RWSIZE;
  Scalar* RWORK = new Scalar[RWORKdim];

  //
  // compute L-infinity norm of matrix
  //
  double ANORM = SymmetricSparse<double>::GetMaxNorm();

  double EPS   = 1e-14;
  char   CMACH = 'E';
  EPS          = dlamch_(&CMACH);
  double TOL   = EPS * ANORM;
  IWSIZE       = 3 * this->d_numRow + 2 * d_Factor->d_NSUPER;

  int ASDEF = 0;

  //       -------------------------------------------------------
  //       BLKLDL ...  numerical factorization.
  //
  //       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX, XLNZ,
  //                   LNZ, DEFBLK, TOL, TMPSIZ, IWSIZE, RWSIZE
  //       Output:     LNZ, NDEF, LBDEF, DEF, IPROW, IPCOL, IFLAG
  //       Work:       TMPVEC(TMPSIZ), IWORK(2*N+2*NSUPER),
  //                   RWORK(RWSIZE)
  //       -------------------------------------------------------

  if constexpr (std::is_same_v<std::decay_t<Scalar>, double>) {
    blkldl_(
        d_Factor->d_NSUPER,
        d_Factor->d_XSUPER,
        SNODE,
        d_Factor->d_XLINDX,
        d_Factor->d_LINDX,
        d_Factor->d_XLNZ,
        d_Factor->d_LNZ,
        d_Factor->d_DEFBLK,
        ASDEF,
        d_Factor->d_NDEF,
        d_Factor->d_LBDEF,
        d_Factor->d_DEF,
        TOL,
        d_Factor->d_IPROW,
        d_Factor->d_IPCOL,
        TMPSIZ,
        TMPVEC,
        IWSIZE,
        IWORK,
        RWSIZE,
        RWORK,
        IFLAG);
  } else {
    exit(-13);
  }

  if (IFLAG != 0) {
    std::cout << "error in call to blkldl in SymmetricSparse::factor" << std::endl;
    std::cout << "BLKLDL IFLAG=" << IFLAG << std::endl;
  }

  //  if (DEFBLK == 0) {
  //    LBDEF=0;
  //    NDEF=0;
  //  }

  int MAXDEF = 10;
  if ((d_Factor->d_NDEF != 0) && (d_Factor->d_NDEF <= MAXDEF)) {
    //
    // compute null space
    //
    d_Factor->d_NS = new Scalar[this->d_numRow * d_Factor->d_NDEF];
    int LDNS       = this->d_numRow;
    if constexpr (std::is_same_v<std::decay_t<Scalar>, double>) {
      blkns_(
          d_Factor->d_NSUPER,
          d_Factor->d_XSUPER,
          d_Factor->d_XLINDX,
          d_Factor->d_LINDX,
          d_Factor->d_XLNZ,
          d_Factor->d_LNZ,
          d_Factor->d_DEFBLK,
          d_Factor->d_NDEF,
          d_Factor->d_LBDEF,
          d_Factor->d_DEF,
          d_Factor->d_IPCOL,
          d_Factor->d_INVP,
          d_Factor->d_NS,
          LDNS,
          RWORK);
    }
  }  // if ((d_Factor->d_NDEF != 0) && (d_Factor->d_NDEF <= MAXDEF))
  delete[] TMPVEC;
  TMPVEC = 0;
  delete[] RWORK;
  RWORK = 0;
}

template <class Scalar>
void
SymmetricSparse<Scalar>::symbolicFactorization(int*& SNODE, int*& IWORK, int& IWSIZE, int& IFLAG)
{
  //
  // solver parameters
  //
  d_Factor->d_DEFBLK = 1;

  //
  int order_opt = 1;
  //
  int NNZA, NADJ, IWMAX, MAXSUP, NTOT;
  int NNZL, NSUB, NLNZ;

  int TMPSIZ, RWSIZE;

  int OPTIONS[8];
  //
  OPTIONS[0] = 0;  // use default values for options
  OPTIONS[1] = 0;
  MAXSUP     = 150;                        // maximum supernode size
  NTOT       = this->d_NNZ;                // total number of nonzeros in matrix
  NADJ       = NTOT - this->d_numRow;      // number of nonzeros in full matrix minus those on diagonal
  NNZA       = NADJ / 2 + this->d_numRow;  // number of nonzeros in lower triangle including diagonal
  IWMAX      = 7 * this->d_numRow + 3;     // dimension for integer working array
  if ((MAXSUP + 2 * this->d_numRow + 1) > IWMAX) IWMAX = MAXSUP + 2 * this->d_numRow + 1;
  if ((3 * this->d_numRow + 2 * MAXSUP) > IWMAX) IWMAX = 3 * this->d_numRow + 2 * MAXSUP;
  IWORK = new int[IWMAX];

  d_Factor->d_LINDX  = new int[NTOT];
  d_Factor->d_PERM   = new int[this->d_numRow];
  d_Factor->d_INVP   = new int[this->d_numRow];
  d_Factor->d_XSUPER = new int[this->d_numRow + 1];
  d_Factor->d_XLINDX = new int[this->d_numRow + 1];
  d_Factor->d_XLNZ   = new int[this->d_numRow + 1];
  d_Factor->d_DEF    = new int[this->d_numRow];
  d_Factor->d_IPROW  = new int[this->d_numRow];
  d_Factor->d_IPCOL  = new int[this->d_numRow];

  int* ADJ    = new int[NTOT];
  int* XADJ   = new int[this->d_numRow + 1];
  int* COLCNT = new int[this->d_numRow];

  //
  // determine adjacency structure of matrix
  // ADJ[XADJ[i]:XADJ[i+1]-1] = dofs adjacent to dof i (but not including
  //                            dof i itself)
  //
  XADJ[0]    = 0;
  int nnzADJ = 0;
  for (int i = 0; i < this->d_numRow; ++i) {
    for (int j = this->d_ROWPTR[i]; j < this->d_ROWPTR[i + 1]; ++j) {
      int col = this->d_COLIDX[j];
      if (col == i) continue;
      ADJ[nnzADJ] = col;
      nnzADJ++;
    }
    XADJ[i + 1] = nnzADJ;
  }  // for (int i = 0; i < d_numRow; ++i)

  //
  // convert COLPTR, ROWIDX, XADJ, and ADJ to Fortran numbering
  //
  for (int i = 0; i <= this->d_numRow; i++) this->d_ROWPTR[i]++;
  for (int i = 0; i < this->d_NNZ; i++) this->d_COLIDX[i]++;
  for (int i = 0; i <= this->d_numRow; i++) XADJ[i]++;
  for (int i = 0; i < nnzADJ; i++) ADJ[i]++;
  //
  NADJ = XADJ[this->d_numRow] - 1;
  //
  for (int i = 0; i <= this->d_numRow; i++) d_Factor->d_XLINDX[i] = XADJ[i];
  for (int i = 0; i < nnzADJ; i++) d_Factor->d_LINDX[i] = ADJ[i];
  //
  // multiple minimum degree ordering
  //
  IWSIZE = 4 * this->d_numRow;
  if (order_opt == 1) {
    ordmmd2_(
        this->d_numRow,
        d_Factor->d_XLINDX,
        d_Factor->d_LINDX,
        d_Factor->d_INVP,
        d_Factor->d_PERM,
        IWSIZE,
        IWORK,
        NSUB,
        IFLAG);
    if (IFLAG != 0) {
      std::cout << "error in call to ordmmd2 in SymmetricSparse::factor" << std::endl;
      std::cout << "ORDMMD2 IFLAG=" << IFLAG << std::endl;
    }
  }
  //
  // Metis ordering
  //
  if (order_opt == 2) {
    // std::cout << " !!! USE Metis Ordering !!! \n\n";
    int numflag = 1;
    metis_nodend(
        &this->d_numRow, d_Factor->d_XLINDX, d_Factor->d_LINDX, &numflag, OPTIONS, d_Factor->d_PERM, d_Factor->d_INVP);
  }

  // ************************************************************************
  // *****   SFINIT ... set up for symbolic factorization               *****
  // ************************************************************************
  //
  // --------
  // Purpose:
  // --------
  //
  // This subroutine computes the storage requirements and sets up
  // preliminary data structures for the symbolic factorization.
  //
  // Note:
  // This version produces the maximal supernode partition
  // (i.e., the one with the fewest possible supernodes).
  //
  // Input: N, NADJ, XADJ, ADJNCY, PERM, INVP, MAXSUP, DEFBLK, IWSIZE
  //
  // Temporary: IWORK
  //
  // Output: COLCNT, NNZL, NSUB, NSUPER, XSUPER, SNODE, IFLAG
  //

  IWSIZE = 7 * this->d_numRow + 3;
  sfinit_(
      this->d_numRow,
      NADJ,
      XADJ,
      ADJ,
      d_Factor->d_PERM,
      d_Factor->d_INVP,
      MAXSUP,
      d_Factor->d_DEFBLK,
      COLCNT,
      NNZL,
      NSUB,
      d_Factor->d_NSUPER,
      d_Factor->d_XSUPER,
      SNODE,
      IWSIZE,
      IWORK,
      IFLAG);
  //---
  if (IFLAG != 0) {
    std::cout << "error in call to sfinit in SymmetricSparse::factor" << std::endl;
    std::cout << "SFINIT IFLAG=" << IFLAG << std::endl;
  }
  //
  // supernodal symbolic factorization
  //
  IWSIZE = d_Factor->d_NSUPER + 2 * this->d_numRow + 1;
  if (NSUB > NTOT) {
    delete[] d_Factor->d_LINDX;
    d_Factor->d_LINDX = new int[NSUB];
  }
  symfct_(
      this->d_numRow,
      NADJ,
      XADJ,
      ADJ,
      d_Factor->d_PERM,
      d_Factor->d_INVP,
      COLCNT,
      d_Factor->d_NSUPER,
      d_Factor->d_XSUPER,
      SNODE,
      NSUB,
      d_Factor->d_XLINDX,
      d_Factor->d_LINDX,
      d_Factor->d_XLNZ,
      IWSIZE,
      IWORK,
      IFLAG);

  delete[] ADJ;
  delete[] XADJ;
  delete[] COLCNT;

  if (IFLAG != 0) {
    std::cout << "error in call to symfct in SymmetricSparse::factor" << std::endl;
    std::cout << "SYMFCT IFLAG=" << IFLAG << std::endl;
  }
}

template <class Scalar>
int
SymmetricSparse<Scalar>::Solve(int NRHS, Scalar RHS[], Scalar SOL[], Scalar TEMP[])
{
  using NonConstScalar = std::decay_t<Scalar>;
  if constexpr (!std::is_same_v<NonConstScalar, double> && !std::is_same_v<NonConstScalar, std::complex<double>>) {
    std::cout << "\n !!! SymmetricSparse<Scalar>::Solve is not implemented for this type !!!\n\n";
    return -1;
  }

  //
  // numerical solution
  //
  int LRHS  = this->d_numRow;
  int LSOL  = this->d_numRow;
  int LTEMP = this->d_numRow;

  if (this->d_numRow == 1) {
    //--- UH (06/2009)
    // This is a fix when n = 1. There was a bug on the Mac.
    // The modification might not be consistent with the general case.
    //--------
    auto invDiag = 1.0 / d_Factor->d_LNZ[0];
    for (int ic = 0; ic < NRHS; ++ic) { SOL[ic] = RHS[ic] * invDiag; }
  } else {
    //
    //****   BLKSLVN ... Block multiple triangular solutions            *****
    //
    // Given the Cholesky factorization of a sparse symmetric matrix,
    // this subroutine performs the triangular solutions of multiple
    // linear system.  It uses output from BLKLDL.
    //
    if constexpr (std::is_same_v<NonConstScalar, double>) {
      blkslvn_(
          d_Factor->d_NSUPER,
          d_Factor->d_XSUPER,
          d_Factor->d_XLINDX,
          d_Factor->d_LINDX,
          d_Factor->d_XLNZ,
          d_Factor->d_LNZ,
          d_Factor->d_DEFBLK,
          d_Factor->d_NDEF,
          d_Factor->d_LBDEF,
          d_Factor->d_DEF,
          d_Factor->d_IPROW,
          d_Factor->d_IPCOL,
          d_Factor->d_PERM,
          d_Factor->d_INVP,
          LRHS,
          NRHS,
          RHS,
          LSOL,
          SOL,
          LTEMP,
          TEMP);
    } else if constexpr (std::is_same_v<NonConstScalar, std::complex<double>>) {
      zblkslvn_(
          d_Factor->d_NSUPER,
          d_Factor->d_XSUPER,
          d_Factor->d_XLINDX,
          d_Factor->d_LINDX,
          d_Factor->d_XLNZ,
          d_Factor->d_LNZ,
          d_Factor->d_DEFBLK,
          d_Factor->d_NDEF,
          d_Factor->d_LBDEF,
          d_Factor->d_DEF,
          d_Factor->d_IPROW,
          d_Factor->d_IPCOL,
          d_Factor->d_PERM,
          d_Factor->d_INVP,
          LRHS,
          NRHS,
          RHS,
          LSOL,
          SOL,
          LTEMP,
          TEMP);
    }
  }

  return 0;
}

template <class Scalar>
int
SymmetricSparse<Scalar>::Solve(int NRHS, Scalar RHS[], Scalar SOL[])
{
  Scalar* TEMP = new Scalar[NRHS * this->d_numRow];
  Solve(NRHS, RHS, SOL, TEMP);
  delete[] TEMP;
  return 0;
}

template <class Scalar>
SymmetricSparse<Scalar>*
SymmetricSparse<Scalar>::SymmetricExtract(int nFree, int freeToGlobal[]) const
{
  int unordered = 0;
  for (int ii = 1; ii < nFree; ++ii) {
    if (freeToGlobal[ii - 1] >= freeToGlobal[ii]) unordered += 1;
  }
  if (unordered > 0) {
    std::cout << "\n !!! SymmetricSparse<Scalar>::SymmetricExtract >> Index Array is not ordered increasingly !!!\n\n";
  }

  std::vector<int> globalToFree(this->d_numRow, -1);
  for (int ii = 0; ii < nFree; ++ii) { globalToFree[freeToGlobal[ii]] = ii; }

  int my_nnz = 0;
  for (int ii = 0; ii < this->d_numRow; ++ii) {
    if (globalToFree[ii] == -1) continue;
    for (int jj = this->d_ROWPTR[ii]; jj < this->d_ROWPTR[ii + 1]; ++jj) {
      if (globalToFree[this->d_COLIDX[jj]] > -1) my_nnz += 1;
    }
  }

  SymmetricSparse<Scalar>* MyK = new SymmetricSparse<Scalar>(nFree, my_nnz);

  int*    my_rowbeg = MyK->GetRowPointers();
  int*    my_colidx = MyK->GetColIndices();
  Scalar* my_data   = MyK->GetData();

  int count    = 0;
  my_nnz       = 0;
  my_rowbeg[0] = 0;
  for (int ii = 0; ii < this->d_numRow; ++ii) {
    if (globalToFree[ii] == -1) continue;
    my_rowbeg[count + 1] = my_rowbeg[count];
    for (int jj = this->d_ROWPTR[ii]; jj < this->d_ROWPTR[ii + 1]; ++jj) {
      if (globalToFree[this->d_COLIDX[jj]] > -1) {
        my_rowbeg[count + 1] += 1;
        my_colidx[my_nnz] = globalToFree[this->d_COLIDX[jj]];
        my_data[my_nnz]   = this->d_ANZ[jj];
        my_nnz += 1;
      }
    }
    count += 1;
  }

  return MyK;
}

template <class Scalar>
SymmetricSparse<Scalar>*
SymmetricSparse<Scalar>::ExtractSchurComplement(std::vector<int>& keepRows, std::vector<int>& elimRows) const
{
  return this->ExtractSchurComplement(keepRows.size(), &keepRows[0], elimRows.size(), &elimRows[0]);
}

template <class Scalar>
SymmetricSparse<Scalar>*
SymmetricSparse<Scalar>::ExtractSchurComplement(int numKeep, int keepRows[], int numElim, int elimRows[]) const
{
  std::vector<int> gTol_keep(this->GetNumRows(), -1);
  for (int ii = 0; ii < numKeep; ++ii) gTol_keep[keepRows[ii]] = ii;

  SymmetricSparse<Scalar>* M_kk = this->SymmetricExtract(numKeep, keepRows);
  SymmetricSparse<Scalar>* MySC = 0;

  //
  // Get the non-zero indices in the block (elim, keep)
  //
  // If the sparsity of the block (keep, keep) changes,
  // it will be located at those indices.
  //

  std::set<int> sc_fill;
  int*          g_rowbeg = this->GetRowPointers();
  int*          g_colidx = this->GetColIndices();
  for (int ii = 0; ii < numElim; ++ii) {
    int myrow = elimRows[ii];

    for (int jj = g_rowbeg[myrow]; jj < g_rowbeg[myrow + 1]; ++jj) {
      int l_mycol = gTol_keep[g_colidx[jj]];
      if (l_mycol > -1) { sc_fill.insert(l_mycol); }
    }  // for (int jj = g_rowbeg[myrow]; jj < g_rowbeg[myrow+1]; ++jj)

  }  // for (int ii = 0; ii < numElim; ++ii)

  int* mkk_rowbeg = M_kk->GetRowPointers();
  int* mkk_colidx = M_kk->GetColIndices();

  //
  // Create the row pointers for the Schur complement
  //
  int* mysc_rowbeg;
  {
    std::vector<int> sc_rowbeg(numKeep + 1, 0);
    for (int ii = 0; ii < numKeep; ++ii) {
      std::set<int> row_fill;
      //---------
      for (int jj = mkk_rowbeg[ii]; jj < mkk_rowbeg[ii + 1]; ++jj) row_fill.insert(mkk_colidx[jj]);
      for (std::set<int>::iterator jt = sc_fill.begin(); jt != sc_fill.end(); ++jt) row_fill.insert(*jt);
      //---------
      sc_rowbeg[ii + 1] = row_fill.size();
    }

    for (int ii = 1; ii <= numKeep; ++ii) { sc_rowbeg[ii] += sc_rowbeg[ii - 1]; }

    MySC = new SymmetricSparse(numKeep, sc_rowbeg[numKeep]);

    //
    // Copy the row pointers
    //
    mysc_rowbeg = MySC->GetRowPointers();
    for (int ii = 0; ii <= numKeep; ++ii) mysc_rowbeg[ii] = sc_rowbeg[ii];
  }

  //
  // Create the column indices
  //
  int* mysc_colidx = MySC->GetColIndices();
  for (int ii = 0; ii < numKeep; ++ii) {
    std::set<int> row_fill;
    //---------
    for (int jj = mkk_rowbeg[ii]; jj < mkk_rowbeg[ii + 1]; ++jj) row_fill.insert(mkk_colidx[jj]);
    for (std::set<int>::iterator jt = sc_fill.begin(); jt != sc_fill.end(); ++jt) row_fill.insert(*jt);
    //---------
    int count = mysc_rowbeg[ii];
    for (std::set<int>::iterator jt = row_fill.begin(); jt != row_fill.end(); ++jt) mysc_colidx[count++] = *jt;
  }

  //
  // Add the entries from M_kk
  //
  Scalar* mkk_values  = M_kk->GetData();
  Scalar* mysc_values = MySC->GetData();

  for (int ii = 0; ii < numKeep; ++ii) {
    for (int jj = mkk_rowbeg[ii]; jj < mkk_rowbeg[ii + 1]; ++jj) {
      int mycol = mkk_colidx[jj];
      for (int kk = mysc_rowbeg[ii]; kk < mysc_rowbeg[ii + 1]; ++kk) {
        if (mysc_colidx[kk] == mycol) {
          mysc_values[kk] += mkk_values[jj];
          break;
        }
      }
    }
  }

  //
  // Delete the first block for memory
  //
  delete M_kk;

  //
  // Subtract the eliminated part Mke * Mee^{-1} * Mek
  //

  SymmetricSparse<Scalar>* M_ee = this->SymmetricExtract(numElim, elimRows);
  M_ee->Factor();

  SparseMatrix<Scalar>* M_ke = this->Extract(numKeep, keepRows, numElim, elimRows);

  int     blockSize = 4;
  Scalar* TmpSol    = new Scalar[2 * blockSize * numKeep];
  for (int iv = 0; iv < numKeep; iv += blockSize) {
    int     numVec = (iv + blockSize > numKeep) ? numKeep - iv : blockSize;
    Scalar* FullProd;
    {
      SparseMatrix<Scalar>* M_ek = this->Extract(numElim, elimRows, numVec, keepRows + iv);
      FullProd                   = Full(*M_ek);
      delete M_ek;
    }
    //----------
    M_ee->Solve(numVec, FullProd, TmpSol, TmpSol + numVec * numElim);
    //----------
    M_ke->Apply(numVec, TmpSol, FullProd);
    //----------
    //
    // Store the product now
    //
    for (std::set<int>::iterator it = sc_fill.begin(); it != sc_fill.end(); ++it) {
      int row = *it;
      for (int ic = 0; ic < numVec; ++ic) {
        int col = keepRows[iv + ic];
        for (int kk = mysc_rowbeg[row]; kk < mysc_rowbeg[row + 1]; ++kk) {
          if (mysc_colidx[kk] == col) {
            mysc_values[kk] -= FullProd[row + ic * numKeep];
            break;
          }
        }
      }
    }  // for (std::set<int>::iterator it = sc_fill.begin(); it != sc_fill.end(); ++it)

    //--------------
    // For Debugging
    //--------------
    double maxNorm = 0.0;
    for (int itmp = 0; itmp < numVec * numKeep; ++itmp)
      maxNorm = (maxNorm > std::abs(FullProd[itmp])) ? maxNorm : std::abs(FullProd[itmp]);
    for (int ir = 0; ir < numKeep; ++ir) {
      if (sc_fill.find(ir) != sc_fill.end()) continue;
      for (int ic = 0; ic < numVec; ++ic) {
        if (FullProd[ir + ic * numKeep] > 1.0e-15 * maxNorm) {
          std::cout << " !!! LARGE ENTRY " << ir << ", " << keepRows[iv + ic];
          std::cout << " = " << FullProd[ir + ic * numKeep] << " " << maxNorm << std::endl;
        }
      }
    }
    //----------
    delete[] FullProd;
  }
  delete[] TmpSol;
  delete M_ee;
  delete M_ke;

  return MySC;
}

template <class Scalar>
void
SymmetricSparse<Scalar>::inpnv(
    int&   n,
    int    colptr[],
    int    rowidx[],
    Scalar values[],
    int    perm[],
    int    invp[],
    int&   nsuper,
    int    xsuper[],
    int    xlindx[],
    int    lindx[],
    int    xlnz[],
    Scalar lnz[],
    int    offset[])
{
  //       ***************************************************
  //       Numerical input into data structure for sparse LDL'
  //       factorization.
  //       ***************************************************
  //
  // --------
  // Purpose:
  // --------
  //
  // This subroutine inputs numerical values of a symmetric matrix
  // into sparse data structures that have been set up for Cholesky
  // factorization.  It is assumed that the input matrix is stored
  // by columns.
  //
  //       --------------------------------------------------------
  //       INPNV ...   input numerical values into data structures.
  //
  //       Input:      N, COLPTR, ROWIDX, VALUES, PERM, INVP,
  //                   NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
  //       Output:     LNZ
  //       Work:       IWORK(N)
  //       --------------------------------------------------------
  //
  // Note (10/17/08 -- UH):
  // This code is simply a C-translation of the Fortran routine inpnv.f
  //

  int ii, lxbeg, fstcol, jsuper, lxend, jlen, irow, lstcol, jcol, oldj, lastl;

  for (ii = 1; ii < xlnz[n]; ii++) lnz[ii - 1] = (Scalar)0;

  lxbeg  = xlindx[0];
  fstcol = xsuper[0];

  for (jsuper = 1; jsuper <= nsuper; jsuper++) {
    // -----------------------------------------------
    // First get offset to facilitate numerical input.
    // -----------------------------------------------
    lxend = xlindx[jsuper];
    jlen  = lxend - lxbeg;

    for (ii = lxbeg; ii < lxend; ii++) {
      irow = lindx[ii - 1];
      jlen--;
      offset[irow - 1] = jlen;
    }

    lstcol = xsuper[jsuper];

    for (jcol = fstcol; jcol < lstcol; jcol++) {
      // -----------------------------------
      // Next input the individual nonzeros.
      // -----------------------------------
      oldj  = perm[jcol - 1];
      lastl = xlnz[jcol] - 1;
      for (ii = colptr[oldj - 1]; ii < colptr[oldj]; ii++) {
        irow = invp[rowidx[ii - 1] - 1];
        if (irow >= fstcol) { lnz[lastl - offset[irow - 1] - 1] = values[ii - 1]; }  // if (irow >= fstcol)
      }  // for (ii=colptr[oldj-1];ii<colptr[oldj];ii++)
    }  // for (jcol=fstcol;jcol<lstcol;jcol++)

    lxbeg  = lxend;
    fstcol = lstcol;

  }  // for (jsuper=1;jsuper<=nsuper;jsuper++)
}
