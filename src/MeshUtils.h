#pragma once

#include "Element.h"
#include "Mesh.h"

#include <string>
#include <vector>

namespace IMSI {

enum class DomainType: char { InputFile, 
                              Bar, 
                              Rectangle, Trapeze, 
                              Brick
                            };

struct DomainParams {
  std::string fileName = "";
  double lowerCorner[3] = {0.0, 0.0, 0.0};
  double upperCorner[3] = {1.0, 1.0, 1.0};
  int numElePerDir[3] = {0, 0, 0};
  DomainType omega = DomainType::Rectangle;
  ElementType cellType = ElementType::Q1;
};

/// \brief Function to generate the finite element mesh
Mesh GenerateMesh(DomainParams const& params, std::vector<double> corners = {});

    void GetMatrixSparsity(const Mesh& grid);

}
