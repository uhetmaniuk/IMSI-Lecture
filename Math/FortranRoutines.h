#ifndef MATHCXX_FORTRAN_ROUTINES_H
#define MATHCXX_FORTRAN_ROUTINES_H


#include <complex>


#ifdef FEC_USE_MPI
#include "mpi.h"
#endif


namespace MathCXX 
{

// LOGICAL as 4 bytes
typedef int LOGICAL;

  #if defined (INTEL_CXML)
    #ifndef FEC_F77PREFIX
      #define FEC_F77PREFIX __stdcall
    #endif
    #ifndef FEC_F77FUNC
      #define FEC_F77FUNC(lcase,UCASE) UCASE
    #endif
    #ifndef FEC_F77CHAR
      #define FEC_F77CHAR char *, unsigned int
    #endif
  #else
    #ifndef FEC_F77PREFIX
      #define FEC_F77PREFIX
    #endif
    #ifndef FEC_F77FUNC
      #define FEC_F77FUNC(lcase,UCASE) lcase ##_  
    #endif
    #ifndef FEC_F77CHAR
      #define FEC_F77CHAR char *
    #endif
  #endif

  #ifdef CHAR_MACRO
    #undef CHAR_MACRO
  #endif

  #if defined (INTEL_CXML)
    #define CHAR_MACRO(char_var) &char_var, 1
  #else
    #define CHAR_MACRO(char_var) &char_var
  #endif

  #ifdef __cplusplus
  extern "C" {
  #endif

  // Double precision BLAS 1 //
  void FEC_F77PREFIX FEC_F77FUNC(daxpy,DAXPY)(int *n, double *a, double x[], int* incx, 
                                              double y[], int* incy);
  double FEC_F77PREFIX FEC_F77FUNC(ddot,DDOT)(int *n, double x[], int* incx, double y[], 
                                              int* incy);
  double FEC_F77PREFIX FEC_F77FUNC(dnrm2,DNRM2)(int *n, double x[], int* incx);
  void FEC_F77PREFIX FEC_F77FUNC(dscal,DSCAL)(int *n, double *a, double x[], int* incx);
  void FEC_F77PREFIX FEC_F77FUNC(dswap,DSWAP)(int *n, double x[], int* incx, double y[], 
                                              int* incy);

  // Double precision BLAS 2 //
  void FEC_F77PREFIX FEC_F77FUNC(dgemv,DGEMV)(FEC_F77CHAR, int* m, int* n, double* alpha, 
                                 double A[], int* lda, double x[], int* incx, 
                                 double* beta, double y[], int* incy);
  
  // Double precision BLAS 3 //
  void FEC_F77PREFIX FEC_F77FUNC(dgemm,DGEMM)(FEC_F77CHAR, FEC_F77CHAR, int* m, int* n, int* k,
                                 double* alpha, double A[], int* lda, double x[], int* incx, 
                                 double* beta, double y[], int* incy);
  
  // Double precision LAPACK //
  void FEC_F77PREFIX FEC_F77FUNC(dgels,DGELS)(FEC_F77CHAR, int *m, int* n, int* nrhs, double A[],
                                 int* lda, double B[], int *ldb, double *work, int *lwork,
                                 int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dgeqrf,DGEQRF)(int *M, int *N, double *A, int *lda, double *tau, 
                       double *work, int *lwork, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dgetrf,DGETRF)(int *M, int *N, double A[], int *lda, 
                                                int ipiv[], int *info); 
  void FEC_F77PREFIX FEC_F77FUNC(dgetrs,DGETRS)(FEC_F77CHAR, int* n, int* nrhs, double A[],
                                 int* lda, int ipiv[], double B[], int *ldb, int *info); 
  void FEC_F77PREFIX FEC_F77FUNC(dgesvd,DGESVD)(FEC_F77CHAR, FEC_F77CHAR, int *M, int *N, 
                                 double *A, int *lda, double *s, double *u, int *ldu, double *vt,
                                 int *ldvt, double *work, int *lwork, int *info);
  double FEC_F77PREFIX FEC_F77FUNC(dlamch,DLAMCH)(FEC_F77CHAR);
  void FEC_F77PREFIX FEC_F77FUNC(dormqr,DORMQR)(FEC_F77CHAR, FEC_F77CHAR, int *M, int *N, int *K, double *A,
                       int *lda, double *tau, double *C, int *ldc, double *work, int *lwork,
                       int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dpotrf,DPOTRF)(FEC_F77CHAR, int *n, double *a, int *lda, 
                                 int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dpotrs,DPOTRS)(FEC_F77CHAR, int *n, int *nrhs, double *a, 
                                 int *lda, double *x, int *ldx, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dpteqr,DPTEQR)(FEC_F77CHAR, int *N, double *D, double *E,
                       double *Z, int *ldz, double *work, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dspev,DSPEV)(FEC_F77CHAR, FEC_F77CHAR, int *N, double *A, double *W, 
                       double *Z, int *ldz, double *work, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dsteqr,DSTEQR)(FEC_F77CHAR, int *N, double *D, double *E, double *Z,
                       int *ldz, double *work, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dstev,DSTEV)(FEC_F77CHAR, int *N, double *D, double *E, double *Z,
                                              int *ldz, double *work, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dsyev,DSYEV)(FEC_F77CHAR, FEC_F77CHAR, int *N, double *A, int *lda,
                       double *W, double *work, int *lwork, int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dsygv,DSYGV)(int *itype, FEC_F77CHAR, FEC_F77CHAR, int *N, double *A,
                       int *lda, double *B, int *ldb, double *W, double *work, int *lwork,
                       int *info);
  void FEC_F77PREFIX FEC_F77FUNC(dtrsm,DTRSM)(FEC_F77CHAR, FEC_F77CHAR, FEC_F77CHAR, 
                       FEC_F77CHAR, int *M, int *N, double *alpha, double *A, int *lda,
                       double *B, int *ldb);

  
  // Complex Double precision LAPACK //
  void FEC_F77PREFIX FEC_F77FUNC(zgesvd,ZGESVD)(FEC_F77CHAR, FEC_F77CHAR, int *M, int *N, 
    std::complex<double> *A, int *lda, double *s, 
    std::complex<double> *u, int *ldu, std::complex<double>  *vt, int *ldvt, 
    std::complex<double> *work, int *lwork, int *info);

  
#if defined (INTEL_CXML)
  int FEC_F77PREFIX FEC_F77FUNC(ilaenv,ILAENV)(int *ispec, char *NAME, unsigned int len_name, char *OPTS, 
                      unsigned int len_opts, int *N1, int *N2, int *N3, int *N4);
#else
  int FEC_F77PREFIX FEC_F77FUNC(ilaenv,ILAENV)(int *ispec, char *NAME, char *OPTS, int *N1, int *N2,
                      int *N3, int *N4, int len_name, int len_opts);
#endif

    
#ifdef FEC_USE_ARPACK
  // Double precision customized ARPACK routines //
  #if defined (INTEL_CXML)
  void FEC_F77PREFIX FEC_F77FUNC(mydsaupd,MYDSAUPD)(int *, char *, unsigned int, int *, char *, 
                       unsigned int, int *, double *, double *, int *, double *, int *, int *,
                       int *, double *, double *, int *, int *, int *);
  void FEC_F77PREFIX FEC_F77FUNC(dsaupd,DSAUPD)(int *, char *, unsigned int, int *, char *, 
  unsigned int, int *, double *, double *, int *, double *, int *, int *,
  int *, double *, double *, int *, int *, int *);
  void FEC_F77PREFIX FEC_F77FUNC(dseupd,DSEUPD)(LOGICAL *rvec, char *HOWMNY, unsigned int len_howny,
                       LOGICAL *select, double *D, double *Z, int *ldz, double *sigma, char *BMAT,
                       unsigned int len_bmat, int *N, char *which, unsigned int len_which, 
                       int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                       int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
                       int *info);
  #else
  void FEC_F77PREFIX FEC_F77FUNC(mydsaupd,MYDSAUPD)(int *, char *, int *, char *, int *, 
                       double *, double *, int *, double *, int *, int *, int *, double *,
                       double *, int *, int *, int *, int, int);
  void FEC_F77PREFIX FEC_F77FUNC(dsaupd,DSAUPD)(int *, char *, int *, char *, int *, 
  double *, double *, int *, double *, int *, int *, int *, double *,
  double *, int *, int *, int *, int, int);
  void FEC_F77PREFIX FEC_F77FUNC(dseupd,DSEUPD)(LOGICAL *rvec, char *HOWMNY, LOGICAL *select, double *D, 
                       double *Z, int *ldz, double *sigma, char *BMAT, int *N, char *which,
                       int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                       int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
                       int *info, int len_howmny, int len_bmat, int len_which);
  #endif

  #ifdef FEC_USE_MPI
  // Double precision customized PARPACK routines //
  #if defined (INTEL_CXML)
  void FEC_F77PREFIX FEC_F77FUNC(mypdsaupd,MYPDSAUPD)(MPI_Comm *, int *, char *, unsigned int, int *, 
                       char *, unsigned int, int *, double *, double *, int *, double *,
                       int *, int *, int *, double *, double *, int *, int *, int *);
  void FEC_F77PREFIX FEC_F77FUNC(pdsaupd,PDSAUPD)(MPI_Comm *, int *, char *, unsigned int, int *, 
                       char *, unsigned int, int *, double *, double *, int *, double *,
                       int *, int *, int *, double *, double *, int *, int *, int *);
  void FEC_F77PREFIX FEC_F77FUNC(pdseupd,PDSEUPD)(MPI_Comm *MyComm, LOGICAL *rvec, char *HOWMNY, 
                       unsigned int len_howmny, LOGICAL *select, double *D, double *Z,
                       int *ldz, double *sigma, char *BMAT, unsigned int len_bmat, int *N,
                       char *which, unsigned int len_which, int *nev, double *tol, double *resid,
                       int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
                       double *workl, int *lworkl, int *info);
  #else
  void FEC_F77PREFIX FEC_F77FUNC(mypdsaupd,MYPDSAUPD)(MPI_Comm *, int *, char *, int *, char *, int *, 
                       double *, double *, int *, double *, int *, int *, int *, double *,
                       double *, int *, int *, int *, int, int);
  void FEC_F77PREFIX FEC_F77FUNC(pdsaupd,PDSAUPD)(MPI_Comm *, int *, char *, int *, char *, int *, 
  double *, double *, int *, double *, int *, int *, int *, double *,
  double *, int *, int *, int *, int, int);
  void FEC_F77PREFIX FEC_F77FUNC(pdseupd,PDSEUPD)(MPI_Comm *MyComm, LOGICAL *rvec, char *HOWMNY, 
                       LOGICAL *select, double *D, double *Z, int *ldz, double *sigma,
                       char *BMAT, int *N, char *which, int *nev, double *tol, double *resid, 
                       int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
                       double *workl, int *lworkl, int *info, int len_howmny,
                       int len_bmat, int len_which);
  #endif

  // End for FEC_USE_MPI
  #endif

// End for FEC_USE_ARPACK
#endif

  #ifdef __cplusplus
  }
  #endif

  //! A class to interface Fortran routines.

  /// To be done
  /// later.
  class FortranRoutines {

    public: 

    //-------------------------//
    // Double precision BLAS 1 //
    //-------------------------//


    //! Constant times a vector plus a vector.

    /// Wrapper around the BLAS routine 'daxpy'
    /// \note See BLAS documentation (http://www.netlib.org/blas/daxpy.f)
    void AXPY(int N, double ALPHA, double *X, int incx, double *Y, int incy) const;

    //! Constant times a vector plus a vector.

    /// Wrapper around the BLAS routine 'daxpy'
    /// \note See BLAS documentation (http://www.netlib.org/blas/daxpy.f)
    void AXPY(int N, double ALPHA, double *X, double *Y) const
         { AXPY(N, ALPHA, X, 1, Y, 1); }

    //! Forms the Euclidian dot product of two vectors.

    /// Wrapper around the BLAS routine 'ddot'
    /// \note See BLAS documentation (http://www.netlib.org/blas/ddot.f)
    double DOT(int N, double *X, int incx, double *Y, int incy) const;

    //! Forms the Euclidian dot product of two vectors.

    /// Wrapper around the BLAS routine 'ddot'
    /// \note See BLAS documentation (http://www.netlib.org/blas/ddot.f)
    double DOT(int N, double *X, double *Y) const
           { return DOT(N, X, 1, Y, 1); }

    //! Returns the Euclidian norm of a vector.

    /// Wrapper around the BLAS routine 'dnrm2'
    /// \note See BLAS documentation (http://www.netlib.org/blas/dnrm2.f)
    double NRM2(int N, double *X, int incx = 1) const;

    //! Scales a vector by a constant.

    /// Wrapper around the BLAS routine 'dscal'
    /// \note See BLAS documentation (http://www.netlib.org/blas/dscal.f)
    void SCAL(int N, double ALPHA, double *X, int incX = 1) const;

    //! Interchanges two vectors.

    /// Wrapper around the BLAS routine 'dswap'
    /// \note See BLAS documentation (http://www.netlib.org/blas/dswap.f)
    void SWAP(int N, double *X, int incx, double *Y, int incy) const;

    
    //-------------------------//
    // Double precision BLAS 2 //
    //-------------------------//

    
    //!  DGEMV  performs matrix-vector products.

    /// Wrapper around the BLAS routine 'dgemv'
    ///  DGEMV  performs one of the matrix-vector operations
    ///
    ///     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    ///
    ///  where alpha and beta are scalars, x and y are vectors and A is an
    ///  m by n matrix.
    /// \note See BLAS documentation (http://www.netlib.org/blas/dgemv.f)
    void GEMV(char TRANS, int M, int N, double ALPHA, double * A, int LDA, 
              double * X, int INCX, double BETA, double * Y, int INCY) const;

    //!  DGEMV  performs matrix-vector products.

    /// Wrapper around the BLAS routine 'dgemv'
    ///  DGEMV  performs one of the matrix-vector operations
    ///
    ///     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    ///
    ///  where alpha and beta are scalars, x and y are vectors and A is an
    ///  m by n matrix.
    /// \note See BLAS documentation (http://www.netlib.org/blas/dgemv.f)
    void GEMV(char TRANS, int M, int N, double ALPHA, double * A, int LDA, 
              double * X, double BETA, double * Y) const
              { GEMV(TRANS, M, N, ALPHA, A, LDA, X, 1, BETA, Y, 1); }

    //-------------------------//
    // Double precision BLAS 3 //
    //-------------------------//

    //!  DGEMM  performs matrix-matrix operations

    /// Wrapper around the BLAS routine 'dgemm'
    ///  DGEMM  performs one of the matrix-matrix operations
    ///
    ///     C := alpha*op( A )*op( B ) + beta*C,
    ///
    ///  where  op( X ) is one of
    ///
    ///     op( X ) = X   or   op( X ) = X',
    ///
    ///  alpha and beta are scalars, and A, B and C are matrices, with op( A )
    ///  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    /// \note See BLAS documentation (http://www.netlib.org/blas/dgemm.f)
    void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
              double ALPHA, double * A, int LDA, double * B,
              int LDB, double BETA, double * C, int LDC) const;

    //! DTRSM  solves a linear system with a triangular matrix.

    /// Wrapper around the BLAS routine 'dtrsm'
    ///  DTRSM  solves one of the matrix equations
    ///
    ///     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
    ///
    ///  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    ///  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
    ///
    ///     op( A ) = A   or   op( A ) = A'.
    ///
    ///  The matrix X is overwritten on B.
    /// \note See LAPACK documentation (http://www.netlib.org/blas/dtrsm.f)
    void TRSM(char SIDE, char UPLO, char TRANS, char DIAG, int m, int n, double alpha,
              double *A, int lda, double *B, int ldb) const;

    //-----------------//
    // LAPACK routines //
    //-----------------//

    //!  DGELS solves overdetermined or underdetermined real linear systems.

    /// Wrapper around the LAPACK routine 'dgels'
    ///  DGELS solves overdetermined or underdetermined real linear systems
    ///  involving an M-by-N matrix A, or its transpose, using a QR or LQ
    ///  factorization of A.  It is assumed that A has full rank.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dgels.f)
    void GELS(char TRANS, int M, int N, int NRHS, double *A, int lda, double *B, int ldb,
              double *work, int lwork, int *info) const;

    //!  DGEQRF computes a QR factorization of a real M-by-N matrix A: A = Q * R.

    /// Wrapper around the LAPACK routine 'dgeqrf'
    ///  DGEQRF computes a QR factorization of a real M-by-N matrix A: A = Q * R.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dgeqrf.f)
    void GEQRF(int M, int N, double *A, int lda, double *tau, double *work, int lwork, 
               int *info) const;

    //!  DGETRF computes an LU factorization of a general matrix.

    /// Wrapper around the LAPACK routine 'dgetrs'
    ///  DGETRF computes an LU factorization of a general M-by-N matrix A
    ///  using partial pivoting with row interchanges.
    ///
    ///  The factorization has the form
    ///     A = P * L * U
    ///  where P is a permutation matrix, L is lower triangular with unit
    ///  diagonal elements (lower trapezoidal if m > n), and U is upper
    ///  triangular (upper trapezoidal if m < n).
    ///
    ///  This is the right-looking Level 3 BLAS version of the algorithm.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dgetrf.f)
    void GETRF(int m, int n, double *A, int lda, int *ipiv, int *info) const;

    //!  DGETRS solves a system of linear equations

    /// Wrapper around the LAPACK routine 'dgetrs'
    ///  DGETRS solves a system of linear equations
    ///     A * X = B  or  A' * X = B
    ///  with a general N-by-N matrix A using the LU factorization computed
    ///  by DGETRF.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dgetrs.f)
    void GETRS(char TRANS, int n, int nrhs, double *A, int lda, int *ipiv, double *B,
               int ldb, int *info) const;

    //!  DGESVD computes the singular value decomposition (SVD) of a real matrix.

    /// Wrapper around the LAPACK routine 'dgesvd'
    ///  DGESVD computes the singular value decomposition (SVD) of a real
    ///  M-by-N matrix A, optionally computing the left and/or right singular
    ///  vectors. The SVD is written
    ///
    ///       A = U * SIGMA * transpose(V)
    ///
    ///  where SIGMA is an M-by-N matrix which is zero except for its
    ///  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
    ///  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
    ///  are the singular values of A; they are real and non-negative, and
    ///  are returned in descending order.  The first min(m,n) columns of
    ///  U and V are the left and right singular vectors of A.
    ///
    ///  Note that the routine returns V**T, not V.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dgesvd.f)
    void GESVD(char JOBU, char JOBVT, int M, int N, double * A, int LDA, double * S, 
              double * U, int LDU, double * VT, int LDVT, double * WORK, int * LWORK, 
              int * INFO) const;

    //!  ILAENV is called from LAPACK routines to choose problem-dependent parameters.

    /// Wrapper around the LAPACK routine 'ilaenv'
    ///  ILAENV is called from the LAPACK routines to choose problem-dependent
    ///  parameters for the local environment.  See ISPEC for a description of
    ///  the parameters.
    ///
    ///  This version provides a set of parameters which should give good,
    ///  but not optimal, performance on many of the currently available
    ///  computers.  Users are encouraged to modify this subroutine to set
    ///  the tuning parameters for their particular machine using the option
    ///  and problem size information in the arguments.
    ///
    ///  This routine will not function correctly if it is converted to all
    ///  lower case.  Converting it to all upper case is allowed.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/util/ilaenv.f)
    int LAENV(int ispec, char *NAME, char *OPTS, int N1, int N2, int N3, int N4, 
              int len_name, int len_opts) const;

    //!  DLAMCH determines double precision machine parameters.

    /// Wrapper around the LAPACK routine 'dlamch'
    ///  DLAMCH determines double precision machine parameters.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/util/dlamch.f)
    void LAMCH(char CMACH, double &T) const;

    //!  DORMQR performs matrix-matrix multiplications for one orthogonal matrix.

    /// Wrapper around the LAPACK routine 'dormqr'
    ///  DORMQR overwrites the general real M-by-N matrix C with
    ///
    ///                  SIDE = 'L'     SIDE = 'R'
    ///  TRANS = 'N':      Q * C          C * Q
    ///  TRANS = 'T':      Q**T * C       C * Q**T
    ///
    ///  where Q is a real orthogonal matrix defined as the product of k
    ///  elementary reflectors
    ///
    ///        Q = H(1) H(2) . . . H(k)
    ///
    ///  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
    ///  if SIDE = 'R'.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dormqr.f)
    void ORMQR(char SIDE, char TRANS, int M, int N, int K, double *A, int lda, double *tau, 
               double *C, int ldc, double *work, int lwork, int *info) const;

    //!  DPOTRF computes the Cholesky factorization of a real symmetric positive definite matrix.

    /// Wrapper around the LAPACK routine 'dpotrf'
    ///  DPOTRF computes the Cholesky factorization of a real symmetric
    ///  positive definite matrix A.
    ///
    ///  The factorization has the form
    ///     A = U**T * U,  if UPLO = 'U', or
    ///     A = L  * L**T,  if UPLO = 'L',
    ///  where U is an upper triangular matrix and L is lower triangular.
    ///
    ///  This is the block version of the algorithm, calling Level 3 BLAS.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dpotrf.f)
    void POTRF(char UPLO, int N, double *A, int LDA, int *INFO) const;

    //!  DPOTRS solves a system with a positive definite matrix using the Cholesky factorization.

    /// Wrapper around the LAPACK routine 'dpotrs'
    ///  DPOTRS solves a system of linear equations A*X = B with a symmetric
    ///  positive definite matrix A using the Cholesky factorization
    ///  A = U**T*U or A = L*L**T computed by DPOTRF.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dpotrs.f)
    void POTRS(char UPLO, int N, int NRHS, double *A, int LDA, double *X, int LDX, int *INFO) const;

    //! DPTEQR computes the eigenpairs of a symmetric positive definite tridiagonal matrix.

    /// DPTEQR computes all eigenvalues and, optionally, eigenvectors of a
    /// symmetric positive definite tridiagonal matrix by first factoring the
    /// matrix using DPTTRF, and then calling DBDSQR to compute the singular
    /// values of the bidiagonal factor.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dpteqr.f)
    void PTEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, double *work, 
               int *info) const;

    //!  DSPEV computes eigenpairs of a real symmetric matrix in packed storage.

    /// Wrapper around the LAPACK routine 'dspev'
    ///  DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
    ///  real symmetric matrix A in packed storage.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dspev.f)
    void SPEV(char JOBZ, char UPLO, int N, double *A, double *W, double *Z, int ldz, 
              double *work, int *info) const;

    //!  DSTEQR computes eigenpairs of a symmetric tridiagonal matrix.

    /// Wrapper around the LAPACK routine 'dsteqr'
    ///  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
    ///  symmetric tridiagonal matrix using the implicit QL or QR method.
    ///  The eigenvectors of a full or band symmetric matrix can also be found
    ///  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
    ///  tridiagonal form.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dsteqr.f)
    void STEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, double *work,
               int *info) const;

    //!  DSTEV computes eigenpairs of a tridiagonal symmetric matrix.
    
    /// Wrapper around the LAPACK routine 'dstev'
    ///  DSTEV computes all eigenvalues and, optionally, eigenvectors of a
    ///  tridiagonal symmetric matrix A.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dstev.f)
    void STEV(char JOBZ, int N, double *D, double *E, double *Z, int ldz, double *work,
              int *info) const;

    //!  DSYEV computes eigenpairs of a real symmetric matrix.

    /// Wrapper around the LAPACK routine 'dsyev'
    ///  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
    ///  real symmetric matrix A.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dsyev.f)
    void SYEV(char JOBZ, char UPLO, int N, double *A, int lda, double *W, double *work,
              int lwork, int *info) const;

    //!  DSYGV computes eigenpairs of a real generalized symmetric-definite eigenproblem.

    /// Wrapper around the LAPACK routine 'dsygv'
    ///  DSYGV computes all the eigenvalues, and optionally, the eigenvectors
    ///  of a real generalized symmetric-definite eigenproblem, of the form
    ///  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    ///  Here A and B are assumed to be symmetric and B is also
    ///  positive definite.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/double/dsygv.f)
    void SYGV(int itype, char JOBZ, char UPLO, int N, double *A, int lda, double *B, int ldb,
              double *W, double *work, int lwork, int *info) const;

    
    //!  ZGESVD computes the singular value decomposition (SVD) of a complex matrix.

    /// Wrapper around the LAPACK routine 'zgesvd'
    ///  ZGESVD computes the singular value decomposition (SVD) of a complex
    /// M-by-N matrix A, optionally computing the left and/or right singular
    ///  vectors. The SVD is written
    ///
    ///       A = U * SIGMA * conjugate-transpose(V)
    ///
    ///  where SIGMA is an M-by-N matrix which is zero except for its
    ///  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
    ///  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
    ///  are the singular values of A; they are real and non-negative, and
    ///  are returned in descending order.  The first min(m,n) columns of
    ///  U and V are the left and right singular vectors of A.
    ///
    ///  Note that the routine returns V**H, not V.
    /// \note See LAPACK documentation (http://www.netlib.org/lapack/patch-3.0/src/zgesvd.f)
    void GESVD(char JOBU, char JOBVT, int M, int N, 
               std::complex<double> * A, int LDA, double * S, 
               std::complex<double> * U, int LDU, std::complex<double> * VT, int LDVT, 
               std::complex<double> * WORK, int * LWORK, 
               int * INFO) const;

    
    //---------------------------------//
    // Double precision ARPACK routines
    //---------------------------------//

    //! ARPACK DSAUPD
    void SAUPD(int *ido, char BMAT, int N, char *which, int nev, double tol, double *resid,
               int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd, double *workl,
               int lworkl, int *info, int verbose) const;

    //! ARPACK DSEUPD
    void SEUPD(LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, double *Z, int ldz,
               double sigma, char BMAT, int N, char *which, int nev, double tol, double *resid,
               int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
               double *workl, int lworkl, int *info) const;

  #ifdef FEC_USE_MPI
    // Double precision PARPACK routines
    //! ARPACK PDSAUPD
    void PSAUPD(MPI_Comm MyComm, int *ido, char BMAT, int N, char *which, int nev, double tol, 
                double *resid, int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
                double *workl, int lworkl, int *info, int verbose) const;
    //! ARPACK PDSEUPD
    void PSEUPD(MPI_Comm MyComm, LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, double *Z, 
                int ldz, double sigma, char BMAT, int N, char *which, int nev, double tol,
                double *resid, int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
                double *workl, int lworkl, int *info) const;

  #endif

  };

}

#endif
