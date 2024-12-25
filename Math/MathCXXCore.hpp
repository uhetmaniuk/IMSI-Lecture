#ifndef MATH_CXX_CORE_HPP
#define MATH_CXX_CORE_HPP


#include <complex>
#include <string>
#include <vector>


namespace MathCXX {
  
  
  static bool d_hasInitializedSRand = false;
  
  
  //! Function to give a random real value.
  /// \remark The random number is in [-1, 1].
  void SetToRandom(double &val);

  //! Function to give a random real value.
  /// \remark The random number is in [-1, 1].
  void SetToRandom(float &val);

  //! Function to give a random integer value.
  void SetToRandom(int &val);

  //! Function to give a random complex value.
  /// \remark Each random component belongs to [-1, 1].
  void SetToRandom(std::complex<double> &val);

  //! Function to give a random complex value.
  /// \remark Each random component belongs to [-1, 1].
  void SetToRandom(std::complex<float> &val);


  //
  // Defines function involving block vectors
  //
  
  
  template<class T> class BlockVector;
  
  
  /// Defines the component-wise cosine of a block vector.
  template <class T>
  const BlockVector<T> cos
  (
   const BlockVector<T> &x
   );
  
  
  /// Overloads the output operator for BlockVector<T>.
  template<class T> std::ostream& operator<< 
  (
   std::ostream & os, 
   const BlockVector<T> &ref
   );
    

  //
  // Defines functions involving dense matrices
  //
  
  
  template<class T> class DenseMatrix;
  
  
  /// Defines the component-wise absolute value of a matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> abs
  (
    const DenseMatrix<T> &b
  );
  
  
  /// Defines the component-wise absolute value of a complex matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> abs
  (
      const DenseMatrix< std::complex<T> > &b
  );
  
  
  /// Computes the condition number of a matrix.
  /// \note Computes the singular values of the matrix.
  template <class T>
  T cond
  (
    const DenseMatrix< T > &a
  );


  /// Computes the condition number of a complex matrix.
  /// \note Computes the singular values of the matrix.
  template <class T>
  T cond
  (
    const DenseMatrix< std::complex<T> > &a
  );


  /// Computes the determinant of a matrix.
  template <class T>
  T det
  (
    const DenseMatrix< T > &a
  );


  /// Defines the component-wise imaginary part of a matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> imag
  (
    const DenseMatrix<T> &b
  );
  
  
  /// Defines the component-wise imaginary part of a complex matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> imag
  (
    const DenseMatrix< std::complex<T> > &b
  );
  
  
  /// Computes the Frobenius norm of a dense matrix.
  template <class T>
  T normFrobenius
  (
    const DenseMatrix< T > &a
  );


  /// Computes the Frobenius norm of a complex dense matrix.
  template <class T>
  T normFrobenius
  (
    const DenseMatrix< std::complex<T> > &a
  );


  /// Computes the norm of a matrix.
  template <class T>
  T norm
  (
    const DenseMatrix< T > &x, const size_t &p
  );


  /// Computes the norm of a complex matrix.
  template <class T>
  T norm
  (
    const DenseMatrix< std::complex<T> > &x, const size_t &p
  );


  /// Defines the multiplication of a scalar and a matrix.
  template <class T>
  const DenseMatrix<T> operator*
  (
   const T &a,
   const DenseMatrix<T> &b
  );
  
  
  /// Defines the multiplication of a scalar and a complex matrix.
  template <class T>
  const DenseMatrix< std::complex<T> > operator*
  (
   const T &a,
   const DenseMatrix< std::complex<T> > &b
  );
  
  
  /// Defines the multiplication of a complex scalar and a matrix.
  template <class T>
  const DenseMatrix< std::complex<T> > operator*
  (
   const std::complex<T> &a,
   const DenseMatrix<T> &b
  );
  
 
  /// Overloads the output operator for DenseMatrix<T>.
  template<class T> std::ostream& operator<< 
  (
    std::ostream & os, 
    const DenseMatrix<T> &ref
  );
  
  
  /// Defines the component-wise real part of a matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> real
  (
    const DenseMatrix<T> &b
  );
  
  
  /// Defines the component-wise real part of a complex matrix (a la Matlab).
  template <class T>
  const DenseMatrix<T> real
  (
    const DenseMatrix< std::complex<T> > &b
  );
  
  
  /// Computes the singular values of a matrix.
  /// \param[in] B  Reference to a dense matrix.
  /// \param[in] sing  Pointer to the array of singular values.
  template <class T>
  int sv
  (
    const DenseMatrix<T> &B
    , T* sing
  );
  

template <>
int sv<double>
(
  const DenseMatrix<double> &b,
  double *sing
);

  
  /// Computes the singular values of a complex matrix.
  /// \param[in] B  Reference to a dense matrix.
  /// \param[in] sing  Pointer to the array of singular values.
  template <class T>
  int sv
  (
    const DenseMatrix< std::complex<T> > &B
    , T* sing
  );
  

  //
  // Defines functions involving tridiagonal matrices
  //
  
  
  template <class T> class TridiagonalMatrix;
  
  
  /// Computes the eigenvalues for a tridiagonal matrix.
  /// \param[in] B  Reference to a tridiagonal matrix.
  /// \param[in] eval  Pointer to the array of eigenvalues.
  /// \param[in] flag  Flag to exploit structure in the matrix (Default = NONE, SYM_UPPER, SYM_LOWER)
  /// \note The length of 'eval' should be equal to the size of the matrix.
  /// \note When the flag is SYM_UPPER, only the values on the diagonal and the upper diagonal are
  /// accessed.
  template <class T>
  int eig
  (
   const TridiagonalMatrix<T> &B
   , T* eval
   , const std::string flag = "NONE"
   );
  
  
  /// Computes the eigenvalues and eigenvectors for a symmetric tridiagonal matrix.
  /// \param[in] B  Reference to a symmetric tridiagonal matrix.
  /// \param[in] eval  Pointer to the array of eigenvalues.
  /// \param[in] Z  Reference to a dense matrix storing the eigenvectors.
  /// \param[in] flag  Flag to exploit structure in the matrix (Default = NONE, SYM_UPPER, SYM_LOWER)
  /// \note The length of 'eval' should be equal to the size of the matrix.
  /// \note The matrix Z should be of same size as the input matrix.
  /// \note When the flag is SYM_UPPER, only the values on the diagonal and the upper diagonal are
  /// accessed.
  template <class T>
  int eig
  (
   const TridiagonalMatrix<T> &B
   , T* eval
   , DenseMatrix<T> &Z
   , const std::string flag = "NONE"
   );
  
  
  /// Overloads the output operator for TridiagonalMatrix<T>.
  template<class T> std::ostream& operator<<
  (
   std::ostream & os,
   const TridiagonalMatrix<T> &ref
   );
  
  
  //
  // Functions related to quadrature rules
  //
  
  
  //! Function to get quadrature points.
  template <class T>
  void getQuadrature
  (
   const std::string &Type,
   const int sdim,
   const int order,
   std::vector<T> &Points,
   std::vector<T> &Weights
   );
 

}


#endif


