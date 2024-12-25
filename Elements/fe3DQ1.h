#ifndef FECODE_FE_3DQ1_H
#define FECODE_FE_3DQ1_H


#include "Elements/LagrangeElement.h"


namespace FECode {
  
  //! A class for three-dimensional Q1 finite element.
  
  /// The class defines the shape functions for the trilinear isoparametric element.
  ///
  /// The nodes are numbered as follows
  ///
  ///    Front face
  ///
  ///    4 --- 3 
  ///
  ///    1 --- 2
  ///
  ///    Back face
  ///
  ///    8 --- 7 
  ///
  ///    5 --- 6
  ///
  /// We recommend a Gaussian quadrature with, at least, 2 points per direction.
  ///
  class fe3DQ1: public virtual LagrangeElement 
  {
    
  private:
    
    std::vector<double> d_JxW;
    unsigned int d_numGauss;
    std::vector<double> d_Phi;
    std::vector<double> d_PhiGrad;
    std::vector<double> d_PointInt;

    static const int sdim = 3;
    static const int nnode = 8;

    // Array for the point coordinates in the quadrature rule.
    // The length of d_lQuad must be three times the length of d_wQuad.
    std::vector<double> d_lQuad;
    
    // Array for the weights in the quadrature rule
    std::vector<double> d_wQuad;
    
    // Don't define these functions
    fe3DQ1(const fe3DQ1 &ref);
    fe3DQ1 & operator[](const fe3DQ1 &ref);
    
  public:
    
    /// Constructor
    fe3DQ1
    (
     const std::vector<double> &nodeCoord, 
     const std::vector<double> &lQ,
     const std::vector<double> &wQ
     );
    
    
    /// Returns the number of integration points (INLINE).
    unsigned int NumIntegrationPoints() const 
    { return d_numGauss; }; 
    
    
    /// Returns a selected integration point.
    void IntegrationPoint
    (
     const size_t iQ, 
     std::vector<double> &vv
     )
    { 
      vv = std::vector<double>(&d_PointInt[sdim*iQ], &d_PointInt[sdim*iQ] + sdim);
    };
    
    
    /// Returns the polynomial degree for the element (INLINE)
    int PolynomialDegree() const 
    { return 1; };
    
    
    //! Gives the value of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sv  Value of the shape function.
    void GetShapeValue
    (
     const size_t iP, const size_t iQ, double &sv
     )
    { sv = d_Phi[iP * d_numGauss + iQ]; };
    
    
    //! Gives the gradient of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sv  Value of the gradient.
    void GetShapeGrad
    (
     const size_t iP
     , const size_t  iQ
     , std::valarray<double> &sg
     );
    
    
    //! Gives the Hessian value of a specified shape function at a specified integration point.
    /// \param[in] iP  Index for the shape function.
    /// \param[in] iQ  Index for the integration point.
    /// \param[out] sh  Value of the Hessian for the shape function.
    void GetShapeHessian
    (
     const size_t iP, const size_t  iQ
     , std::vector<double> &sh
     );
    
    
    /// Returns the volume of the element.
    double Volume();
    
    //! Returns the value of the weighted Jacobian at a specified quadrature point.
    /// \param[in] iQ  Index for the specified quadrature point.
    /// \return Value of the weighted Jacobian.
    double WeightedJacobian
    (
     unsigned int iQ
     )
    { return d_JxW[iQ]; };
    
  };
  
}

#endif
