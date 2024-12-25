/* 
For INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
#if defined (INTEL_CXML)
        unsigned int one=1;
#endif
*/


#include "FortranRoutines.h"


using namespace MathCXX;


// Double precision BLAS 1 //


void FortranRoutines::AXPY(int N, double ALPHA, double *X, int incx, double *Y, int incy) const 
{
  FEC_F77FUNC(daxpy,DAXPY)(&N, &ALPHA, X, &incx, Y, &incy);
}

double FortranRoutines::DOT(int N, double *X, int incx, double *Y, int incy) const 
{
  return FEC_F77FUNC(ddot,DDOT)(&N, X, &incx, Y, &incy);
}

double FortranRoutines::NRM2(int N, double *X, int incx) const 
{
  return FEC_F77FUNC(dnrm2,DNRM2)(&N, X, &incx);
}

void FortranRoutines::SCAL(int N, double ALPHA, double *X, int incX) const 
{
  FEC_F77FUNC(dscal,DSCAL)(&N, &ALPHA, X, &incX);
}

void FortranRoutines::SWAP(int N, double *X, int incx, double *Y, int incy) const 
{
  FEC_F77FUNC(dswap,DSWAP)(&N, X, &incx, Y, &incy);
}


// Double precision BLAS 2 //


void FortranRoutines::GEMV(char TRANS, int M, int N, double ALPHA, double * A, int LDA, 
                           double * X, int INCX, double BETA, double * Y, int INCY) const 
{
  FEC_F77FUNC(dgemv,DGEMV)(CHAR_MACRO(TRANS),&M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}


// Double precision BLAS 3 //


void FortranRoutines::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
                           double ALPHA, double * A, int LDA, double * X, int INCX,
                           double BETA, double * Y, int INCY) const 
{
  FEC_F77FUNC(dgemm,DGEMM)(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, 
                           &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}


// Double precision LAPACK //


void FortranRoutines::GELS(char TRANS, int M, int N, int NRHS, double *A, int lda,
                           double *B, int ldb, double *work, int lwork, int *info) const 
{
  FEC_F77FUNC(dgels,DGELS)(CHAR_MACRO(TRANS), &M, &N, &NRHS, A, &lda, B, &ldb, 
                           work, &lwork, info);
}

void FortranRoutines::GEQRF(int M, int N, double *A, int lda, double *tau, double *work,
                            int lwork, int *info) const 
{
  FEC_F77FUNC(dgeqrf,DGEQRF)(&M, &N, A, &lda, tau, work, &lwork, info);
}

void FortranRoutines::GETRF(int m, int n, double *A, int lda, int *ipiv, int *info) const 
{
  FEC_F77FUNC(dgetrf,DGETRF)(&m, &n, A, &lda, ipiv, info);
}

void FortranRoutines::GETRS(char TRANS, int n, int nrhs, double *A, int lda, int *ipiv, 
                            double *B, int ldb, int *info) const 
{
  FEC_F77FUNC(dgetrs,DGETRS)(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, ipiv, B, &ldb, info);
}

void FortranRoutines::GESVD(char JOBU, char JOBVT, int M, int N, double * A, int LDA, 
                            double * S, double * U, int LDU, double * VT, int LDVT, 
                            double * WORK, int * LWORK, int * INFO) const 
{
  FEC_F77FUNC(dgesvd,DGESVD)(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &M, &N, A, &LDA,
                             S, U, &LDU, VT, &LDVT, WORK, LWORK, INFO);
}

void FortranRoutines::GESVD(char JOBU, char JOBVT, int M, int N, 
                            std::complex<double> * A, int LDA, 
                            double * S, 
                            std::complex<double> * U, int LDU, std::complex<double> * VT, int LDVT, 
                            std::complex<double> * WORK, int * LWORK, int * INFO) const 
{
  FEC_F77FUNC(zgesvd,ZGESVD)(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &M, &N, A, &LDA,
  S, U, &LDU, VT, &LDVT, WORK, LWORK, INFO);
}

int FortranRoutines::LAENV(int ispec, char *NAME, char *OPTS, int N1, int N2, int N3,
                           int N4, int len_name, int len_opts) const 
{
#if defined (INTEL_CXML)
  return FEC_F77FUNC(ilaenv,ILAENV)(&ispec, NAME, len_name, OPTS, len_opts, &N1, &N2, &N3, &N4);
#else
  return FEC_F77FUNC(ilaenv,ILAENV)(&ispec, NAME, OPTS, &N1, &N2, &N3, &N4, len_name, len_opts);
#endif
}

void FortranRoutines::LAMCH(char CMACH, double &T) const 
{
  T = FEC_F77FUNC(dlamch,DLAMCH)(CHAR_MACRO(CMACH));
}

void FortranRoutines::ORMQR(char SIDE, char TRANS, int M, int N, int K, double *A, int lda, 
                            double *tau, double *C, int ldc, double *work, int lwork,
                            int *info) const 
{
  FEC_F77FUNC(dormqr,DORMQR)(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &K, A, &lda, tau, 
                          C, &ldc, work, &lwork, info);
}

void FortranRoutines::POTRF(char UPLO, int N, double *A, int LDA, int *INFO) const 
{
  FEC_F77FUNC(dpotrf,DPOTRF)(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}

void FortranRoutines::POTRS(char UPLO, int N, int NRHS, double *A, int LDA, 
                            double *X, int LDX, int *INFO) const 
{
  FEC_F77FUNC(dpotrs,DPOTRS)(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}

void FortranRoutines::PTEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, 
                            double *work, int *info) const 
{
  FEC_F77FUNC(dpteqr,DPTEQR)(CHAR_MACRO(COMPZ), &N, D, E, Z, &ldz, work, info);
}

void FortranRoutines::SPEV(char JOBZ, char UPLO, int N, double *A, double *W, double *Z,
                           int ldz, double *work, int *info) const 
{
  FEC_F77FUNC(dspev,DSPEV)(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, W, Z, &ldz, work, info);
}

void FortranRoutines::STEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, 
                            double *work, int *info) const 
{
  FEC_F77FUNC(dsteqr,DSTEQR)(CHAR_MACRO(COMPZ), &N, D, E, Z, &ldz, work, info);
}

void FortranRoutines::STEV(char JOBZ, int N, double *D, double *E, double *Z, int ldz, 
                           double *work, int *info) const
{
  FEC_F77FUNC(dstev,DSTEV)(CHAR_MACRO(JOBZ), &N, D, E, Z, &ldz, work, info);
}

void FortranRoutines::SYEV(char JOBZ, char UPLO, int N, double *A, int lda, double *W,
                           double *work, int lwork, int *info) const 
{
  FEC_F77FUNC(dsyev,DSYEV)(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &lda, W, work, &lwork, info);
}

void FortranRoutines::SYGV(int itype, char JOBZ, char UPLO, int N, double *A, int lda, 
                           double *B, int ldb, double *W, double *work, int lwork, 
                           int *info) const 
{
  FEC_F77FUNC(dsygv,DSYGV)(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &lda, B, &ldb,
           W, work, &lwork, info);
}


void FortranRoutines::TRSM(char SIDE, char UPLO, char TRANS, char DIAG, int m, int n, 
                           double alpha, double *A, int lda, double *B, int ldb) const 
{
  FEC_F77FUNC(dtrsm,DTRSM)(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), 
                           CHAR_MACRO(DIAG), &m, &n, &alpha, A, &lda, B, &ldb);
}


// Double precision ARPACK routines


void FortranRoutines::SAUPD(int *ido, char BMAT, int N, char *which, int nev, double tol, 
                            double *resid, int ncv, double *V, int ldv, int *iparam, 
                            int *ipntr, double *workd, double *workl, int lworkl, int *info,
                            int verbose) const 
{
#ifdef FEC_USE_ARPACK
#if defined (INTEL_CXML)
//  FEC_F77FUNC(mydsaupd,MYDSAUPD)(ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, V, &ldv,
//           iparam, ipntr, workd, workl, &lworkl, info, &verbose);
  FEC_F77FUNC(dsaupd,DSAUPD)(ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose);
#else
//  FEC_F77FUNC(mydsaupd,MYDSAUPD)(ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
//           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
  FEC_F77FUNC(dsaupd,DSAUPD)(ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
#endif
#endif
}

void FortranRoutines::SEUPD(LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, 
                            double *Z, int ldz, double sigma, char BMAT, int N, 
                            char *which, int nev, double tol, double *resid, int ncv, double *V,
                            int ldv, int *iparam, int *ipntr, double *workd, double *workl,
                            int lworkl, int *info) const 
{
#ifdef FEC_USE_ARPACK
#if defined (INTEL_CXML)
  FEC_F77FUNC(dseupd,DSEUPD)(&rvec, &HOWMNY, 1, select, D, Z, &ldz, &sigma, &BMAT, 1, &N,
           which, 2, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl,
           info);
#else
  FEC_F77FUNC(dseupd,DSEUPD)(&rvec, &HOWMNY, select, D, Z, &ldz, &sigma, &BMAT, &N, which, 
           &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, 
           1, 1, 2);
#endif
#endif
}

#ifdef FEC_USE_MPI

// Double precision PARPACK routines

void FortranRoutines::PSAUPD(MPI_Comm MyComm, int *ido, char BMAT, int N, char *which, int nev, 
                             double tol, double *resid, int ncv, double *V, int ldv, int *iparam,
                             int *ipntr, double *workd, double *workl, int lworkl, int *info, 
                             int verbose) const 
{
#ifdef FEC_USE_ARPACK
#if defined (INTEL_CXML)
//  FEC_F77FUNC(mypdsaupd,MYPDSAUPD)(&MyComm, ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, 
//           V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, &verbose);
  FEC_F77FUNC(pdsaupd,PDSAUPD)(&MyComm, ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, 
           V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, &verbose);
#else
//  FEC_F77FUNC(mypdsaupd,MYPDSAUPD)(&MyComm, ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
//           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
  FEC_F77FUNC(pdsaupd,PDSAUPD)(&MyComm, ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
#endif
#endif
}

void FortranRoutines::PSEUPD(MPI_Comm MyComm, LOGICAL rvec, char HOWMNY, LOGICAL *select, 
                             double *D, double *Z, int ldz, double sigma, char BMAT, int N,
                             char *which, int nev, double tol, double *resid, int ncv, double *V,
                             int ldv, int *iparam, int *ipntr, double *workd, double *workl,
                             int lworkl, int *info) const 
{
#ifdef FEC_USE_ARPACK
#if defined (INTEL_CXML)
  FEC_F77FUNC(pdseupd,PDSEUPD)(&MyComm, &rvec, &HOWMNY, 1, select, D, Z, &ldz, &sigma, &BMAT, 1, &N,
           which, 2, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
#else
  FEC_F77FUNC(pdseupd,PDSEUPD)(&MyComm, &rvec, &HOWMNY, select, D, Z, &ldz, &sigma, &BMAT, &N,
           which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, 
           1, 1, 2);
#endif
#endif
}

#endif

