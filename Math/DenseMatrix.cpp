#include "DenseMatrix.hpp"
#include "FortranRoutines.h"

#include <iostream>



/////////////////////////////////////
//
// Definition of specialized member functions
//
////////////////////////////////////


namespace MathCXX {
    
  //----------------------------------------------------------
  
  template <>
  int DenseMatrix< double >::Apply
  (
   const BlockVector<double> &X,
   BlockVector<double> &Y
  )
  {
    
    std::cout << " ... DenseMatrix<double>::Apply --> Specialization \n";
    
    if ((Y.NumRows() != d_numRow) || (X.NumRows() != d_numCol)
        || (Y.NumVectors() < X.NumVectors()))
      return -1;
    
    double *xval = X.Values();
    double *yval = Y.Values();
    size_t xcol = X.NumVectors();
    size_t xrow = X.NumRows();
    
    MathCXX::FortranRoutines call;
    call.GEMM('N', 'N', d_numRow, xcol, xrow, 1.0, 
              d_val_p, d_leadDim, xval, xrow,
              0.0, yval, d_numRow);
    
    return 0;
    
  }
  
  //----------------------------------------------------------
  
}


