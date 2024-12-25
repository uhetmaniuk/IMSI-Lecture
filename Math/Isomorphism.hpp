#ifndef MATHCXX_ISOMORPHISM_HPP
#define MATHCXX_ISOMORPHISM_HPP


#include "Morphism.hpp"


namespace MathCXX 
{
  
  template <class Scalar> class BlockVector;
  
  //! A templated abstract class to store an isomorphism operator.
  
  /// This class defines an interface for creating a mathematical isomorphism.
  template<class Scalar>
  class Isomorphism: public virtual Morphism<Scalar> {
    
  public:
    
    //! Destructor.
    virtual ~Isomorphism
    (
    ) 
    { };
    
    //! Solve a block of systems with the isomorphism.
    //! \param[in] X: Input block vector for the operator.
    //! \param[out] Y: Output block vector after solving the system.
    virtual int Solve
    (
     const BlockVector<Scalar> &X, 
     BlockVector<Scalar> &Y
     ) = 0;
    
  };
  
}

#endif

