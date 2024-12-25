#include "Geometry/Mesh.h"

#include <cassert>
#include <ostream>
#include <iostream>

#include "Geometry/GeometricElement.h"

namespace msfem {

Mesh::~Mesh
()
{
    for (auto &cellPtr : d_cellList) {
        delete cellPtr;
    }
}

  std::ostream & operator<<
  (
   std::ostream & os, const msfem::Mesh &ref
  ) 
  {
    
    os << "//\n// NODES\n//\n";
    
      for (size_t iN = 0; iN < ref.d_nodeData.size(); iN += d_dim) {
      os.precision(8);
      os.setf(std::ios::scientific, std::ios::floatfield);
      for (int jD = 0; jD < d_dim; ++jD)
        os << ref.d_nodeData[iN + jD] << " ";
      os << std::endl;
    }
    
    os << "//\n// MESH CONNECTIVITY\n//\n";
    
    for (size_t iE = 0; iE < ref.d_cellList.size(); ++iE) {
        auto list = ref.d_cellList[iE]->NodesList();
        for (auto node : list) {
            os << list[jN] << " ";
        }
      os << std::endl;
    }
    
    return os;
    
  }
  
} // namespace msfem
