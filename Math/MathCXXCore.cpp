#include "BlockVector.hpp"
#include "DenseMatrix.hpp"
#include "FortranRoutines.h"
#include "MathCXXCore.hpp"
#include "Matrix.hpp"
#include "TridiagonalMatrix.hpp"


#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/time.h>


namespace MathCXX 
{
  
//--------------------------------------------------------  
  
  void SetToRandom(double &val)
  {
    double rnd = (double) std::rand() / RAND_MAX;
    val = (double)(-1.0 + 2.0*rnd);
  }
  
  
  void SetToRandom(float &val)
  {
    float rnd = (float) std::rand() / RAND_MAX;
    val = (float)(-1.0 + 2.0*rnd);
  }
  
  
  void SetToRandom(int &val)
  {
    val = std::rand();
  }
  
  
  void SetToRandom(std::complex<double> &val)
  {
    double rnd1 = (double) std::rand() / RAND_MAX;
    double rnd2 = (double) std::rand() / RAND_MAX;
    val = std::complex<double>(-1.0 + 2.0*rnd1, -1.0 + 2.0*rnd2);
  }
  
  
  void SetToRandom(std::complex<float> &val)
  {
    float rnd1 = (float) std::rand() / RAND_MAX;
    float rnd2 = (float) std::rand() / RAND_MAX;
    val = std::complex<float>(-1.0 + 2.0*rnd1, -1.0 + 2.0*rnd2);
  }
  
//--------------------------------------------------------
  
  template <class T>
  const BlockVector<T> cos
  (
   const BlockVector<T> &x
   )
  {
    BlockVector<T> res(x.NumRows(), x.NumVectors());
    const T *xval_p = x.Values();
    size_t xlen = x.NumRows() * x.NumVectors();
    T *rval_p = res.Values();
    for (size_t ii = 0; ii < xlen; ++ii)
      rval_p[ii] = std::cos(xval_p[ii]);
    return res;
  }
  
  template const BlockVector<double> cos<double>(const BlockVector<double> &x);
  template const BlockVector< std::complex<double> > cos< std::complex<double> >
          (const BlockVector< std::complex<double> > &x);
  
//--------------------------------------------------------
  
  template<class T>
  std::ostream& operator<<
  (
   std::ostream & os 
   , const BlockVector<T> &ref
   )
  {
    
    size_t rows = ref.NumRows();
    size_t vecs = ref.NumVectors();
    
    os << " --- Block Vector of dimension ";
    os << rows << "x" << vecs;
    os << " --- " << std::endl;
    for (size_t iR = 0; iR < rows; ++iR) {
      for (size_t jV = 0; jV < vecs; ++jV)
        os << ref(iR, jV) << " ";
      os << std::endl;
    }
    
    return os;
    
  }

  template std::ostream& operator<< <double> 
  (std::ostream & os, const BlockVector<double> &ref);

  template std::ostream& operator<< < std::complex<double> > 
  (std::ostream & os, const BlockVector< std::complex<double> > &ref);
  
//--------------------------------------------------------
  
  template <class T>
  const DenseMatrix<T> abs
  (
   const DenseMatrix<T> &b
   )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const T *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = std::abs(bval_p[pos]);
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> abs<double>(const DenseMatrix<double> &b);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> abs
  (
   const DenseMatrix< std::complex<T> > &b
   )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const std::complex<T> *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = std::abs(bval_p[pos]);
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> abs<double>
      (const DenseMatrix< std::complex<double> > &b);
  
//--------------------------------------------------------

  template <class T>
  T cond
  (
   const DenseMatrix< T > &a
   )
  {
    
    T res = 0;
    
    size_t arow = a.NumRows();
    size_t acol = a.NumColumns();
    size_t amin = (arow > acol) ? acol : arow;
    std::vector<T> sing(amin);
    
    int info = sv(a, &sing[0]);
    if (info == 0)
    {
      res = sing[0] / sing[amin-1];
    }
    
    return res;
    
  }

  template double cond<double>
      (const DenseMatrix<double> &a);
  
//--------------------------------------------------------
  
  template <class T>
  T cond
  (
    const DenseMatrix< std::complex<T> > &a
   )
  {
    
    T res = 0;
    
    size_t arow = a.NumRows();
    size_t acol = a.NumColumns();
    size_t amin = (arow > acol) ? acol : arow;
    std::vector<T> sing(amin);
    
    int info = sv(a, &sing[0]);
    if (info == 0)
    {
      res = sing[0] / sing[amin-1];
    }
    
    return res;
    
  }

  template double cond<double>
      (const DenseMatrix< std::complex<double> > &a);
  
//--------------------------------------------------------
  
  template <class T>
  T det
  (
   const DenseMatrix< T > &a
   )
  {
    
    if (a.NumRows() != a.NumColumns())
      return (T) 0;
    
    T dval = (T) 0;
    switch (a.NumRows())
    {
      case 1:
        dval = a(0,0);
        break;
      case 2:
        dval = a(0,0) * a(1,1) - a(0,1) * a(1,0);
        break;
      case 3:
        {
          dval += a(0,0) * (a(1,1) * a(2,2) - a(1,2) * a(2,1));
          dval -= a(1,0) * (a(0,1) * a(2,2) - a(0,2) * a(2,1));
          dval += a(2,0) * (a(0,1) * a(1,2) - a(1,1) * a(0,2));
        }
        break;
      default:
        std::cerr << "\n The determinant is not implemented for this dimension\n\n";
        break;
    }
    
    return dval;
    
  }

  template double det<double>
      (const DenseMatrix<double> &a);
  
  template std::complex<double> det< std::complex<double> >
      (const DenseMatrix< std::complex<double> > &a);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> imag
  (
    const DenseMatrix<T> &b
  )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const T *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = (T) 0;
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> imag<double>(const DenseMatrix<double> &b);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> imag
  (
    const DenseMatrix< std::complex<T> > &b
  )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const std::complex<T> *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = std::imag(bval_p[pos]);
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> imag<double>
      (const DenseMatrix< std::complex<double> > &b);
  
//--------------------------------------------------------

  template <class T>
  T normFrobenius
  (
   const DenseMatrix< T > &a
   )
  {
    
    T res = (T) 0;
    size_t arow = a.NumRows();
    size_t acol = a.NumColumns();
    for (size_t jj = 0; jj < acol; ++jj)
    {
      for (size_t ii = 0; ii < arow; ++ii)
      {
        res += std::abs(a(ii, jj)) * std::abs(a(ii,jj));
      }
    }
    
    return std::sqrt(res);
    
  }

  template double normFrobenius<double>
      (const DenseMatrix<double> &a);
  
//--------------------------------------------------------

  template <class T>
  T normFrobenius
  (
   const DenseMatrix< std::complex<T> > &a
   )
  {
    
    T res = (T) 0;
    size_t arow = a.NumRows();
    size_t acol = a.NumColumns();
    for (size_t jj = 0; jj < acol; ++jj)
    {
      for (size_t ii = 0; ii < arow; ++ii)
      {
        res += std::abs(a(ii, jj)) * std::abs(a(ii,jj));
      }
    }
    
    return std::sqrt(res);
    
  }

  template double normFrobenius<double>
      (const DenseMatrix< std::complex<double> > &a);
  
//--------------------------------------------------------

  template <class T>
  T norm
  (
   const DenseMatrix<T> &x,
   const size_t &p
   )
  {
    
    size_t xcol = x.NumColumns();
    size_t xrow = x.NumRows();
    T xnorm = (T) 0;
    if (p == 0)
    {
      //--- Compute the max norm
      for (size_t ii = 0; ii < xrow; ++ii)
      {
        T mynorm = (T) 0;
        for (size_t jj = 0; jj < xcol; ++jj)
          mynorm += std::abs(x(ii,jj));
        xnorm = (xnorm > mynorm) ? xnorm : mynorm;
      }
    }
    else if (p == 1)
    {
      //--- Compute the 1-norm
      for (size_t jj = 0; jj < xcol; ++jj)
      {
        T mynorm = (T) 0;
        for (size_t ii = 0; ii < xrow; ++ii)
          mynorm += std::abs(x(ii,jj));
        xnorm = (xnorm > mynorm) ? xnorm : mynorm;
      }
    }
    else
    {
      std::cerr << "\n !!!";
      std::cerr << " The matrix norm is not implemented for this index (Use p = 0).";
      std::cerr << " !!!\n";
      xnorm = norm(x, 0);
    }
    
    return xnorm;
    
  }

  template double norm<double>
      (const DenseMatrix<double> &x, const size_t &p);
  
//--------------------------------------------------------

  template <class T>
  T norm
  (
   const DenseMatrix< std::complex<T> > &x,
   const size_t &p
   )
  {
    
    size_t xcol = x.NumColumns();
    size_t xrow = x.NumRows();
    T xnorm = (T) 0;
    if (p == 0)
    {
      //--- Compute the max norm
      for (size_t ii = 0; ii < xrow; ++ii)
      {
        T mynorm = (T) 0;
        for (size_t jj = 0; jj < xcol; ++jj)
          mynorm += std::abs(x(ii,jj));
        xnorm = (xnorm > mynorm) ? xnorm : mynorm;
      }
    }
    else if (p == 1)
    {
      //--- Compute the 1-norm
      for (size_t jj = 0; jj < xcol; ++jj)
      {
        T mynorm = (T) 0;
        for (size_t ii = 0; ii < xrow; ++ii)
          mynorm += std::abs(x(ii,jj));
        xnorm = (xnorm > mynorm) ? xnorm : mynorm;
      }
    }
    else
    {
      std::cerr << "\n !!!";
      std::cerr << " The matrix norm is not implemented for this index (Use p = 0).";
      std::cerr << " !!!\n";
      xnorm = norm(x, 0);
    }
    
    return xnorm;
    
  }

  template double norm<double>
      (const DenseMatrix< std::complex<double> > &x, const size_t &p);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> operator*
  (
   const T &a,
   const DenseMatrix<T> &b
   )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const T *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = a * bval_p[pos];
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> operator* <double>
      (const double &a, const DenseMatrix<double> &b);
  
  template const DenseMatrix< std::complex<double> > operator* < std::complex<double> >
      (const std::complex<double> &a, const DenseMatrix< std::complex<double> > &b);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix< std::complex<T> > operator*
  (
   const T &a,
   const DenseMatrix< std::complex<T> > &b
   )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix< std::complex<T> > res(brow, bcol, blead);
    const std::complex<T> *bval_p = b.Values();
    std::complex<T> *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = ((std::complex<T>) a) * bval_p[pos];
      }
    }
    
    return res;
    
  }

  template const DenseMatrix< std::complex<double> > operator* <double>
      (const double &a, const DenseMatrix< std::complex<double> > &b);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix< std::complex<T> > operator*
  (
   const std::complex<T> &a,
   const DenseMatrix< T > &b
   )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix< std::complex<T> > res(brow, bcol, blead);
    const T *bval_p = b.Values();
    std::complex<T> *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = a * ((std::complex<T>) bval_p[pos]);
      }
    }
    
    return res;
    
  }

  template const DenseMatrix< std::complex<double> > operator* <double>
      (const std::complex<double> &a, const DenseMatrix<double> &b);
  
//--------------------------------------------------------


  // Overloads the output operator for MathCXX::DenseMatrix<T>.
  template<class T>
  std::ostream& operator<<
  (
   std::ostream & os, 
   const DenseMatrix<T> &ref
   )
  {
    
    size_t rows = ref.NumRows();
    size_t cols = ref.NumColumns();
    
    os << " --- Dense Matrix of dimension ";
    os << rows << "x" << cols;
    os << " --- " << std::endl;
    for (size_t iR = 0; iR < rows; ++iR) {
      for (size_t jV = 0; jV < cols; ++jV)
        os << ref(iR, jV) << " ";
      os << std::endl;
    }
    
    return os;
    
  }

  template std::ostream& operator<< <double> 
      (std::ostream & os, const DenseMatrix<double> &ref);

  template std::ostream& operator<< < std::complex<double> > 
      (std::ostream & os, const DenseMatrix< std::complex<double> > &ref);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> real
  (
    const DenseMatrix<T> &b
  )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const T *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = bval_p[pos];
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> real<double>(const DenseMatrix<double> &b);
  
//--------------------------------------------------------

  template <class T>
  const DenseMatrix<T> real
  (
    const DenseMatrix< std::complex<T> > &b
  )
  {
    
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    DenseMatrix<T> res(brow, bcol, blead);
    const std::complex<T> *bval_p = b.Values();
    T *rval_p = res.Values();
    size_t pos;
    for (size_t jj = 0; jj < bcol; ++jj)
    {
      for (size_t ii = 0; ii < brow; ++ii)
      {
        pos = ii + jj * blead;
        rval_p[pos] = std::real(bval_p[pos]);
      }
    }
    
    return res;
    
  }

  template const DenseMatrix<double> real<double>
      (const DenseMatrix< std::complex<double> > &b);

  
//--------------------------------------------------------


  template <class T>
  int sv
  (
    const DenseMatrix<T> &b,
    T *sing
  )
  {
    
    std::cout << "\n !!! Computing singular values is not implemented for this type !!!\n\n";
    return -1;
    
  }

    
//--------------------------------------------------------
  

  template <>
  int sv<double>
  (
    const DenseMatrix<double> &b,
    double *sing
  )
  {
    
    FortranRoutines call;
    
    int lwork;
    size_t brow = b.NumRows();
    size_t bcol = b.NumColumns();
    size_t blead = b.LeadingDimension();
    
    if (brow < bcol)
    {
      lwork = 3 * brow + bcol;
      lwork = (lwork < 5*brow) ? 5*brow : lwork;
    }
    else
    {
      lwork = 3 * bcol + brow;
      lwork = (lwork < 5*bcol) ? 5*bcol : lwork;
    }
    int minSize = lwork + brow * bcol;
    
    double *work_p = new double[minSize];
    double *bval_p = b.Values();
    
    int info = 0;
    //--- Compute only the singular values
    for (size_t ic = 0; ic < bcol; ++ic)
      memcpy(work_p + ic*brow, bval_p + ic*blead, brow*sizeof(double));
    
    call.GESVD('N', 'N', brow, bcol, work_p, brow,
               sing, (double*) 0, 1, (double*) 0, 1,
               work_p + brow * bcol, &lwork, &info
              );
    
    delete[] work_p;
    
    return info;
    
  }
  

//--------------------------------------------------------
    
  
  template <class T>
  int sv
  (
   const DenseMatrix< std::complex<T> > &B,
   T *sing
   )
  {
    
    std::cout << "\n !!! Computing singular values is not implemented for this type !!!\n\n";
    return -1;
    
  }
  
  
//--------------------------------------------------------

  
  // Overloads the output operator for MathCXX::TridiagonalMatrix<T>.
  template<class T>
  std::ostream& operator<<
  (
   std::ostream & os, 
   const TridiagonalMatrix<T> &ref
   )
  {
    
    size_t rows = ref.NumRows();
    size_t cols = ref.NumColumns();
    
    os << " --- Tridiagonal Matrix of dimension ";
    os << rows << "x" << cols;
    os << " --- " << std::endl;
    for (size_t iR = 0; iR < rows; ++iR) 
    {
      for (size_t jV = 0; jV < cols; ++jV)
        os << ref(iR, jV) << " ";
      os << std::endl;
    }
    
    return os;
    
  }
  
  
  template std::ostream& operator<< <double> 
  (std::ostream & os, const TridiagonalMatrix<double> &ref);
  
  
  template std::ostream& operator<< < std::complex<double> > 
  (std::ostream & os, const TridiagonalMatrix< std::complex<double> > &ref);
  
  
  //--------------------------------------------------------
  
  
  template <class T>
  int eig
  (
   const TridiagonalMatrix<T> &B
   , T *eval
   , const std::string flag
   )
  {
    
    std::cout << "\n !!! Computing eigenvalues is not implemented for this scalar type !!!\n\n";
    return -1;
    
  }
  
  
  //--------------------------------------------------------
  
  
  template <>
  int eig
  (
   const TridiagonalMatrix< double > &B
   , double *eval_p
   , const std::string flag
   )
  {
    
    FortranRoutines call;
    
    int lwork;
    size_t brow = B.NumRows();
    size_t bcol = B.NumColumns();
    
    if ((brow <= 0) || (brow != bcol))
      return -1;
    
    double *Z_p = 0;
    int ldz = 1;
    double *work_p = 0;
    
    int info = 0;
    if (flag == "SYM_UPPER")
    {
      //--- Case for a symmetric tridiagonal matrix with values stored in upper diagonal
      double *B2_p = new double[brow - 1];
      double *bval_p = B.Values();
      for (size_t ii = 0; ii < brow; ++ii)
      {
        eval_p[ii] = bval_p[ii];
        if (ii + 1 < brow)
          B2_p[ii] = bval_p[ii + brow];
      }
      call.STEV('N', brow, eval_p, B2_p, Z_p, ldz, work_p, &info);
      delete[] B2_p;
    }
    else if (flag == "SYM_LOWER")
    {
      //--- Case for a symmetric tridiagonal matrix with values stored in lower diagonal
      double *B2_p = new double[brow - 1];
      double *bval_p = B.Values();
      for (size_t ii = 0; ii < brow; ++ii)
      {
        eval_p[ii] = bval_p[ii];
        if (ii + 1 < brow)
          B2_p[ii] = bval_p[2*brow - 1 + ii];
      }
      call.STEV('N', brow, eval_p, B2_p, Z_p, ldz, work_p, &info);
      delete[] B2_p;
    }
    else
    {
      std::cout << "\n !!! Computing eigenvalues is not implemented for this matrix type !!!\n\n";
      info = -1;
    }
    
    return info;
    
  }
  
  
  //--------------------------------------------------------
  
  
  template <class T>
  int eig
  (
   const TridiagonalMatrix<T> &B
   , T *eval
   , DenseMatrix<T> &Z
   , const std::string flag
   )
  {
    
    std::cout << "\n !!! Computing eigenpairs is not implemented for this scalar type !!!\n\n";
    return -1;
    
  }
  
  
  //--------------------------------------------------------
  
  
  template <>
  int eig
  (
   const TridiagonalMatrix< double > &B
   , double *eval_p
   , DenseMatrix< double > &Z
   , const std::string flag
   )
  {
    
    size_t brow = B.NumRows();
    size_t bcol = B.NumColumns();
    if ((brow <= 0) || (brow != bcol))
      return -1;
    
    if ((Z.NumRows() != brow) || (Z.NumColumns() != brow))
      return -1;
    
    double *Z_p = Z.Values();
    int ldz = Z.LeadingDimension();
    int info = 0;
    
    if (flag == "SYM_UPPER")
    {
      double *work_p = new double[3*brow - 2];
      double *B2_p = work_p + 2*brow - 1;
      double *bval_p = B.Values();
      for (size_t ii = 0; ii < brow; ++ii)
      {
        eval_p[ii] = bval_p[ii];
        if (ii + 1 < brow)
          B2_p[ii] = bval_p[ii + brow];
      }      
      //---------
      FortranRoutines call;
      call.STEV('V', brow, eval_p, B2_p, Z_p, ldz, work_p, &info);
      //---------
      delete[] work_p;
    }
    else if (flag == "SYM_LOWER")
    {
      double *work_p = new double[3*brow - 2];
      double *B2_p = work_p + 2*brow - 1;      
      double *bval_p = B.Values();
      for (size_t ii = 0; ii < brow; ++ii)
      {
        eval_p[ii] = bval_p[ii];
        if (ii + 1 < brow)
          B2_p[ii] = bval_p[2*brow - 1 + ii];
      }
      //---------
      FortranRoutines call;
      call.STEV('V', brow, eval_p, B2_p, Z_p, ldz, work_p, &info);
      //---------
      delete[] work_p;
    }
    else
    {
      std::cout << "\n !!! Computing eigenpairs is not implemented for this matrix type !!!\n\n";
      info = -1;
    }
    
    return info;
    
  }
  
  
  //--------------------------------------------------------

  
  template <class T>
  void getQuadrature
  (
   const std::string &Type,
   const int sdim,
   const int order,
   std::vector<T> &Points,
   std::vector<T> &Weights
  )
  {
    
    std::cout << "\n !!! MathCXX::getQuadrature >> ";
    std::cout << " Not implemented for this scalar type !!!\n\n";
    assert(0 > 1);
    
  }
  
  
  template <>
  void getQuadrature<double>
  (
   const std::string &Type,
   const int sdim,
   const int order,
   std::vector<double> &Points,
   std::vector<double> &Weights
   ) {

	  if ((sdim == 1) && (Type == "Gauss")) {

		  Points.resize(order);
		  Weights.resize(order);

		  switch (order) {
			  case 1: {
				  Points[0] = 0.0;
				  Weights[0] = 2.0;
				  break;
			  }
			  case 2: {
				  Points[0] = -1.0 / sqrt(3.0);
				  Points[1] = -Points[0];
				  //----
				  Weights[0] = 1.0;
				  Weights[1] = 1.0;
				  break;
			  }
			  case 3: {
				  Points[0] = -sqrt(0.6);
				  Points[1] = 0.0;
				  Points[2] = -Points[0];
				  Weights[0] = 5.0 / 9.0;
				  Weights[1] = 8.0 / 9.0;
				  Weights[2] = Weights[0];
				  break;
			  }
			  case 4: {
				  Points[0] = -0.861136311594053;
				  Points[1] = -0.339981043584856;
				  Points[2] = -Points[1];
				  Points[3] = -Points[0];
				  Weights[0] = 0.347854845137454;
				  Weights[1] = 0.652145154862546;
				  Weights[2] = Weights[1];
				  Weights[3] = Weights[0];
				  break;
			  }
			  case 5: {
				  Points[0] = -0.906179845938664;
				  Points[1] = -0.538469310105683;
				  Points[2] = 0.0;
				  Points[3] = -Points[1];
				  Points[4] = -Points[0];
				  Weights[0] = 0.236926885056189;
				  Weights[1] = 0.478628670499366;
				  Weights[2] = 0.568888888888889;
				  Weights[3] = Weights[1];
				  Weights[4] = Weights[0];
				  break;
			  }
			  case 6: {
				  Points[0] = -0.932469514203152;
				  Points[1] = -0.661209386466265;
				  Points[2] = -0.238619186083197;
				  Points[3] = -Points[2];
				  Points[4] = -Points[1];
				  Points[5] = -Points[0];
				  Weights[0] = 0.171324492379170;
				  Weights[1] = 0.360761573048139;
				  Weights[2] = 0.467913934572691;
				  Weights[3] = Weights[2];
				  Weights[4] = Weights[1];
				  Weights[5] = Weights[0];
				  break;
			  }
			  case 7: {
				  Points[0] = -0.949107912342759;
				  Points[1] = -0.741531185599394;
				  Points[2] = -0.405845151377397;
				  Points[3] = 0.0;
				  Points[4] = -Points[2];
				  Points[5] = -Points[1];
				  Points[6] = -Points[0];
				  Weights[0] = 0.129484966168870;
				  Weights[1] = 0.279705391489277;
				  Weights[2] = 0.381830050505119;
				  Weights[3] = 0.417959183673469;
				  Weights[4] = Weights[2];
				  Weights[5] = Weights[1];
				  Weights[6] = Weights[0];
				  break;
			  }
			  case 8: {
				  Points[0] = -0.960289856497536;
				  Points[1] = -0.796666477413627;
				  Points[2] = -0.525532409916329;
				  Points[3] = -0.183434642495650;
				  Points[4] = -Points[3];
				  Points[5] = -Points[2];
				  Points[6] = -Points[1];
				  Points[7] = -Points[0];
				  Weights[0] = 0.101228536290376;
				  Weights[1] = 0.222381034453374;
				  Weights[2] = 0.313706645877887;
				  Weights[3] = 0.362683783378362;
				  Weights[4] = Weights[3];
				  Weights[5] = Weights[2];
				  Weights[6] = Weights[1];
				  Weights[7] = Weights[0];
				  break;
			  }
			  case 9: {
				  Points[0] = -0.968160239507626;
				  Points[1] = -0.836031107326636;
				  Points[2] = -0.613371432700590;
				  Points[3] = -0.324253423403809;
				  Points[4] = 0.0;
				  Points[5] = -Points[3];
				  Points[6] = -Points[2];
				  Points[7] = -Points[1];
				  Points[8] = -Points[0];
				  Weights[0] = 0.081274388361574;
				  Weights[1] = 0.180648160694857;
				  Weights[2] = 0.260610696402935;
				  Weights[3] = 0.312347077040003;
				  Weights[4] = 0.330239355001260;
				  Weights[5] = Weights[3];
				  Weights[6] = Weights[2];
				  Weights[7] = Weights[1];
				  Weights[8] = Weights[0];
				  break;
			  }
			  case 10: {
				  Points[0] = -0.973906528517172;
				  Points[1] = -0.865063366688985;
				  Points[2] = -0.679409568299024;
				  Points[3] = -0.433395394129247;
				  Points[4] = -0.148874338981631;
				  Points[5] = -Points[4];
				  Points[6] = -Points[3];
				  Points[7] = -Points[2];
				  Points[8] = -Points[1];
				  Points[9] = -Points[0];
				  Weights[0] = 0.066671344308688;
				  Weights[1] = 0.149451349150581;
				  Weights[2] = 0.219086362515982;
				  Weights[3] = 0.269266719309996;
				  Weights[4] = 0.295524224714753;
				  Weights[5] = Weights[4];
				  Weights[6] = Weights[3];
				  Weights[7] = Weights[2];
				  Weights[8] = Weights[1];
				  Weights[9] = Weights[0];
				  break;
			  }
			  default: {
				  TridiagonalMatrix<double> BB(order, order);
				  BB = 0.0;
				  for (size_t ii = 1; ii < order; ++ii) {
					  double uu = ((double) ii) / sqrt(4.0 * ii * ii - 1.0);
					  BB(ii - 1, ii) = uu;
				  }
				  DenseMatrix<double> XX(order, order);
				  eig(BB, &Points[0], XX, "SYM_UPPER");
				  for (size_t ii = 0; ii < order; ++ii) {
					  double xrow = XX(0, ii);
					  Weights[ii] = 2.0 * xrow * xrow;
				  }
				  if (order % 2 == 1)
					  Points[order / 2] = 0.0;
				  break;
			  }
		  } // switch(order)

	  } else if ((sdim == 1) && (Type == "GaussLobatto")) {

		  int myOrder = (order > 1) ? order : 2;

		  Points.resize(myOrder);
		  Weights.resize(myOrder);

		  switch (myOrder) {
			  case 2: {
				  Points[0] = -1.0;
				  Points[1] = 1.0;
				  Weights[0] = 1.0;
				  Weights[1] = 1.0;
				  break;
			  }
			  case 3: {
				  Points[0] = -1.0;
				  Points[1] = 0.0;
				  Points[2] = 1.0;
				  Weights[0] = 1.0 / 3.0;
				  Weights[1] = 4.0 * Weights[0];
				  Weights[2] = Weights[0];
				  break;
			  }
			  default: {
				  double *p_old = &Weights[0];
				  double maxDiff = 1.0, eps = 2.5e-16;
				  for (int ii = 0; ii < myOrder; ++ii) {
					  Points[ii] = -std::cos(M_PI * ii / (order - 1.0));
					  p_old[ii] = 2.0;
				  }
				  DenseMatrix<double> PP(myOrder, myOrder);
				  double *pval_p = PP.Values();
				  PP = 1;
				  while (maxDiff > eps) {
					  for (int ii = 0; ii < myOrder; ++ii) {
						  p_old[ii] = Points[ii];
						  pval_p[ii + order] = Points[ii];
					  }
					  for (int kk = 2; kk < myOrder; ++kk) {
						  double ookk = 1.0 / kk;
						  for (int ii = 0; ii < myOrder; ++ii) {
							  PP(ii, kk) = ((2.0 * kk - 1.0) * Points[ii] *
							                PP(ii, kk - 1)
							                - (kk - 1.0) * PP(ii, kk - 2)) *
							               ookk;
						  }
					  }
					  maxDiff = 0.0;
					  for (int ii = 0; ii < myOrder; ++ii) {
						  double pin = PP(ii, myOrder - 1);
						  double tmp =
								  (Points[ii] * pin - PP(ii, myOrder - 2)) /
								  (myOrder * pin);
						  Points[ii] = p_old[ii] - tmp;
						  tmp = std::abs(tmp);
						  maxDiff = (maxDiff > tmp) ? maxDiff : tmp;
					  }
				  } // while (maxDiff > eps)
				  for (int ii = 0; ii < myOrder; ++ii) {
					  double pin = PP(ii, myOrder - 1);
					  Weights[ii] =
							  2.0 / (myOrder * (myOrder - 1.0) * pin * pin);
				  }
				  break;
			  }
		  }

	  }

	  if (sdim == 1)
		  return;

	  if ((Type == "Gauss") || (Type == "GaussLobatto")) {

		  std::vector<double> p1d, w1d;
		  getQuadrature<double>(Type, 1, order, p1d, w1d);

		  if (sdim == 2) {
			  size_t myLen = w1d.size() * w1d.size();
			  Points.resize(2 * myLen);
			  Weights.resize(myLen);
			  size_t counter = 0;
			  for (size_t jq = 0; jq < w1d.size(); ++jq) {
				  for (size_t iq = 0; iq < w1d.size(); ++iq) {
					  Weights[counter] = w1d[iq] * w1d[jq];
					  Points[2 * counter] = p1d[iq];
					  Points[2 * counter + 1] = p1d[jq];
					  counter += 1;
				  }
			  }
		  } else if (sdim == 3) {
			  size_t myLen = w1d.size() * w1d.size() * w1d.size();
			  Points.resize(3 * myLen);
			  Weights.resize(myLen);
			  size_t counter = 0;
			  for (size_t kq = 0; kq < w1d.size(); ++kq) {
				  for (size_t jq = 0; jq < w1d.size(); ++jq) {
					  for (size_t iq = 0; iq < w1d.size(); ++iq) {
						  Weights[counter] = w1d[iq] * w1d[jq] * w1d[kq];
						  Points[3 * counter] = p1d[iq];
						  Points[3 * counter + 1] = p1d[jq];
						  Points[3 * counter + 2] = p1d[kq];
						  counter += 1;
					  }
				  }
			  }
		  } else {
			  std::cout << "\n !!! MathCXX::getQuadrature >> ";
			  std::cout << " Dimension " << sdim
			            << " is not implemented !!!\n\n";
			  assert(0 > 1);
		  }
		  return;

	  }

	  if (Type == "StrangTri") {

	  	if (sdim != 2) {
		    std::cerr << "\n !!! MathCXX::getQuadrature >> ";
		    std::cerr << " 'StrangTri' is only available in 2D !!!\n\n";
		    assert(0 > 1);
	  	}

	  	/* --- 2D Quadrature rule for triangle {(0,0), (1,0), (0,1)} */
	  	switch (order) {
	  		default:
	  		case 1:
	  			Points.assign(2, 1.0/3.0);
	  			Weights.assign(1, 1.0);
	  			break;
	  		case 3:
	  			Points.assign(6, 1.0/6.0);
	  			Points[0] = 2.0/3.0;
	  			Points[3] = Points[0];
			    Weights.assign(3, 1.0/3.0);
	  			break;
	  	}

	  }

	  if (Type == "Keast") {

		  if (sdim != 3) {
			  std::cerr << "\n !!! MathCXX::getQuadrature >> ";
			  std::cerr << " 'Keast' is only available in 3D !!!\n\n";
			  assert(0 > 1);
		  }

		  /* --- 3D Quadrature rule for tetrahedron
		   * {(0,0,0), (1,0,0), (0,1,0), (0,0,1)}
		   */
		  switch (order) {
			  default:
			  case 1:
				  Points.assign(3, 1.0/4.0);
				  Weights.assign(1, 1.0);
				  break;
			  case 4:
				  Points.assign(12, 0.1381966011250105);
				  Points[0] = 0.5854101966249685;
				  Points[4] = Points[0];
				  Points[8] = Points[0];
				  Weights.assign(4, 1.0/4.0);
				  break;
		  }

	  }

  }
  
//--------------------------------------------------------  
  
  
}

