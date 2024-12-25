#ifndef FECODE_RECTANGLE_DOMAIN_H
#define FECODE_RECTANGLE_DOMAIN_H


#include "Geometry/GeometricDomain.h"


namespace FECode {


  class Communicator;
  class ParameterList;
  

  //! A concrete class generating a rectangular domain in two dimensions.

  /// The Rectangle_Domain class defines a geometric domain [a_x, b_x] x [a_y, b_y].
  class Rectangle_Domain: public virtual GeometricDomain 
  {

    protected:

      int d_numEleX;
      int d_numEleY;
    
      double d_ax, d_bx;
      double d_ay, d_by;

      //--- Protected Functions
      void MakeRectangularMesh
      (
        const ParameterList &ParamPb
      );

      // Partitions the domain into subdomains.
      // \note There is one subdomain per processor.
      // \note The routine destroys copies of the global mesh and
      // keeps only local component of the mesh.
      // \note When the number of processors is 4^n, the number of elements 
      // in X-direction is a*2^n, and the number of elements in Y-direction is b*2^n, 
      // the partition is made manually.
      // Else the routine calls METIS to partition the domain.
      // When METIS is not linked to the code, a linear distribution of elements
      // is done.
      void PartitionDomain();

    public:

      /////////////
      // ACCESSOR 
      /////////////

      ////////////
      // CREATOR
      ////////////

      /// Constructor
      /// \param[in] Comm Reference to a Communicator object.
      /// \param[in] ParamPb Reference to a ParameterList object.
      Rectangle_Domain
      (
        const Communicator &Comm,
        const ParameterList &ParamPb
      );

      /// Destructor.
      ~Rectangle_Domain();

      ///////////////
      // MANIPULATOR
      ///////////////

  };

}

#endif

