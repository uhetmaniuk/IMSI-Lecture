#ifndef MATHCXX_MORPHISM_HPP
#define MATHCXX_MORPHISM_HPP


#ifdef FEC_USE_TRILINOS
class Epetra_Operator;
#endif


#include <iostream>


namespace MathCXX 
{

  template <class Scalar> class BlockVector;

  //! A templated abstract class to store a morphism (or operator).
  
  /// The Morphism class enables the use of operators on vectors.
  template <class Scalar>
  class Morphism 
  {
 
    public:
    
      //! Destructor.
      virtual ~Morphism() { };

      //! Apply the operator to a block of vectors.
      //! \param[in] X Block vectors on which the morphism is applied.
      //! \param[in] Y Block vectors storing the result.
      virtual int Apply
      (
        const BlockVector<Scalar> &X, BlockVector<Scalar> &Y
      ) = 0;

#ifdef FEC_USE_TRILINOS
      //! Returns a pointer to an Epetra_Operator object when Trilinos is used.
      //! \return Pointer to Epetra_Operator object
      virtual Epetra_Operator* const getEpetraOperator() 
      { 
        std::cout << "\n !!! Morphism::getEpetraOperator >> Not Implemented !!!\n\n";
        return 0; 
      }; 
#endif

  };

}

#endif

