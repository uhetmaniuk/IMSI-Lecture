#ifndef FECODE_LAGRANGE_ELEMENT_H
#define FECODE_LAGRANGE_ELEMENT_H


#include <valarray>
#include <vector>


namespace FECode {
  
  //! A pure virtual class for Lagrange finite elements.
  
  //! The LagrangeElement class defines the common interface
  //! for nodal-based volumic Lagrange elements.
  class LagrangeElement {
    
  protected:
    
    // Computes the inverse of a matrix n x n stored in J.
    // The input matrix should be stored column-wise.
    // The inverse is then stored in J.
    // Jac returns the determinant of J.
    void inverseJacobian(int n, double *J, double &Jac);
    
  public:
    
    //! Destructor
    virtual ~LagrangeElement() { };
    
    /// Returns the number of integration points.
    //
    // Note that some quadrature rules are specified as one-dimensional.
    // Then the finite element creates the tensor product.
    //
    virtual unsigned int NumIntegrationPoints() const = 0; 
    
    /// Returns the integration points.
    /// \param[in] iQ  Index for the integration point.
    /// \param[in] vv  Vector with coordinates for the integration point.
    virtual void IntegrationPoint
    (
     const size_t iQ
     , std::vector<double> &vv
     ) = 0;
    
    /// Returns the polynomial degree for the element.
    virtual int PolynomialDegree() const { return 0; };
    
    //! Gives the value of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sv  Value of the shape function.
    virtual void GetShapeValue
    (
     const size_t iP, const size_t  iQ
     , double &sv
     ) = 0;
    
    //! Gives the gradient of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sv  Value of the gradient.
    virtual void GetShapeGrad
    (
     const size_t iP
     , const size_t  iQ
     , std::valarray<double> &sg
     ) = 0;


    //! Gives the Hessian value of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sh  Value of the Hessian for the shape function stored in an array.
    virtual void GetShapeHessian
    (
     const size_t iP, const size_t  iQ
     , std::vector<double> &sh
     ) = 0;
    
    
    //! Returns the volume of the element.
    virtual double Volume() = 0;
    
    
    //! Gives the value of the weighted Jacobian at a specified integration point.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] wJ  Value of the weighted Jacobian.
    virtual double WeightedJacobian(unsigned int iQ) = 0;
    
    
    //! Returns the values of gradient of shape functions at a nodal point (INLINE).
    /// \param[in] iP  Index for the nodal point.
    virtual std::vector<double> NodalGrad(unsigned int iP)
    { std::vector<double> res; return res; }
    
    
  };
  
}

#endif
