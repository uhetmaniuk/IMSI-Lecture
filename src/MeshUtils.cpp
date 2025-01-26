#include "MeshUtils.h"

#include <iostream>

#include "FunctionExamples.h"
#include "ScaledLaplacian.h"

namespace IMSI {

    template<unsigned int dim>
    Mesh GenerateUniformTensor(DomainParams const &params) {

        static_assert((dim > 0) && (dim <= 3));

        int const nX = params.numElePerDir[0];
        int const nY = (dim < 2) ? 0 : params.numElePerDir[1];
        int const nZ = (dim < 3) ? 0 : params.numElePerDir[2];
        int const numCells = (dim == 1) ? nX : ((dim == 2) ? nX * nY: nX * nY * nZ);
        auto const Lx = params.upperCorner[0] - params.lowerCorner[0];
        auto const Ly = (dim < 2) ? 0 : params.upperCorner[1] - params.lowerCorner[1];
        auto const Lz = (dim < 3) ? 0 : params.upperCorner[2] - params.lowerCorner[2];
        auto const hX = Lx / nX;
        auto const hY = (dim < 2) ? 0 : Ly / nY;
        auto const hZ = (dim < 3) ? 0 : Lz / nZ;

        // Define global cells
        std::vector<std::vector<int> > cList;
        cList.reserve(numCells);
        switch (params.cellType) {
            case ElementType::Q2: {
                if constexpr (dim == 1) {
                    std::vector<int> list(3);
                    for (int iX = 0; iX < nX; ++iX) {
                        int left = 2 * iX;
                        list[0] = left;
                        list[1] = left + 2;
                        list[2] = left + 1;
                        cList.push_back(list);
                    }
                } else if constexpr (dim == 2) {
                    std::vector<int> list(9);
                    for (int iY = 0; iY < nY; ++iY) {
                        for (int iX = 0; iX < nX; ++iX) {
                            int bottomLeft = 2 * iX + (2 * nX + 1) * 2 * iY;
                            list[0] = bottomLeft;
                            list[1] = bottomLeft + 2;
                            list[2] = bottomLeft + 4 * nX + 4;
                            list[3] = bottomLeft + 4 * nX + 2;
                            list[4] = bottomLeft + 1;
                            list[5] = bottomLeft + 2 * nX + 3;
                            list[6] = bottomLeft + 4 * nX + 3;
                            list[7] = bottomLeft + 2 * nX + 1;
                            list[8] = bottomLeft + 2 * nX + 2;
                            cList.push_back(list);
                        }
                    }
                } else if constexpr (dim == 3) {
                    std::vector<int> list(27);
                    for (int iZ = 0; iZ < nZ; ++iZ) {
                        for (int iY = 0; iY < nY; ++iY) {
                            for (int iX = 0; iX < nX; ++iX) {
                                int bottomLeft = 2 * iX + (2 * nX + 1) * 2 * iY + (2 * nX + 1) * (2 * nY + 1) * 2 * iZ;
                                int middleLeft = bottomLeft + (2 * nX + 1) * (2 * nY + 1);
                                int topLeft = middleLeft + (2 * nX + 1) * (2 * nY + 1);
                                //----------------
                                list[0] = bottomLeft;
                                list[1] = bottomLeft + 2;
                                list[2] = bottomLeft + 4 * nX + 4;
                                list[3] = bottomLeft + 4 * nX + 2;
                                //----------------
                                list[4] = topLeft;
                                list[5] = topLeft + 2;
                                list[6] = topLeft + 4 * nX + 4;
                                list[7] = topLeft + 4 * nX + 2;
                                //----------------
                                list[8] = bottomLeft + 1;
                                list[9] = bottomLeft + 2 * nX + 3;
                                list[10] = bottomLeft + 4 * nX + 3;
                                list[11] = bottomLeft + 2 * nX + 1;
                                //----------------
                                list[12] = topLeft + 1;
                                list[13] = topLeft + 2 * nX + 3;
                                list[14] = topLeft + 4 * nX + 3;
                                list[15] = topLeft + 2 * nX + 1;
                                //----------------
                                list[16] = middleLeft;
                                list[17] = middleLeft + 2;
                                list[18] = middleLeft + 4 * nX + 4;
                                list[19] = middleLeft + 4 * nX + 2;
                                //----------------
                                list[20] = bottomLeft + 2 * nX + 2;
                                //----------------
                                list[21] = topLeft + 2 * nX + 2;
                                //----------------
                                list[22] = middleLeft + 1;
                                //----------------
                                list[23] = middleLeft + 4 * nX + 3;
                                //----------------
                                list[24] = middleLeft + 2 * nX + 1;
                                //----------------
                                list[25] = middleLeft + 2 * nX + 3;
                                //----------------
                                list[26] = middleLeft + 2 * nX + 2;
                                //----------------
                                cList.push_back(list);
                            }
                        }
                    }
                }
                break;
            }
            default: {
                if constexpr (dim == 1) {
                    for (int iX = 0; iX < nX; ++iX) {
                        int bottomLeft = iX;
                        std::vector<int> list(2);
                        list[0] = bottomLeft;
                        list[1] = bottomLeft + 1;
                        cList.push_back(list);
                    }
                } else if constexpr (dim == 2) {
                    for (int iY = 0; iY < nY; ++iY) {
                        for (int iX = 0; iX < nX; ++iX) {
                            int bottomLeft = iX + (nX + 1) * iY;
                            std::vector<int> list(4);
                            list[0] = bottomLeft;
                            list[1] = bottomLeft + 1;
                            list[2] = bottomLeft + nX + 2;
                            list[3] = bottomLeft + nX + 1;
                            cList.push_back(list);
                        }
                    }
                } else if constexpr (dim == 3) {
                    for (int iZ = 0; iZ < nZ; ++iZ) {
                        for (int iY = 0; iY < nY; ++iY) {
                            for (int iX = 0; iX < nX; ++iX) {
                                int bottomLeft = iX + (nX + 1) * iY + (nX + 1) * (nY + 1) * iZ;
                                std::vector<int> list(8);
                                list[0] = bottomLeft;
                                list[1] = bottomLeft + 1;
                                list[2] = bottomLeft + nX + 2;
                                list[3] = bottomLeft + nX + 1;
                                bottomLeft += (nX + 1) * (nY + 1);
                                list[4] = bottomLeft;
                                list[5] = bottomLeft + 1;
                                list[6] = bottomLeft + nX + 2;
                                list[7] = bottomLeft + nX + 1;
                                cList.push_back(list);
                            }
                        }
                    }
                }
                break;
            }
        }

        std::vector<ElementType> cType(numCells, params.cellType);

        // Define global points
        std::vector<std::array<double, 3> > vList;
        std::array<double, 3> coord{0, 0, 0};
        std::vector<int> boundaryNode;
        int nNode = 2;
        switch (params.cellType) {
            case ElementType::Q2: {
                auto const twoNX = 2 * nX;
                auto const twoNY = 2 * nY;
                auto const twoNZ = 2 * nZ;
                vList.reserve((twoNX + 1) * (twoNY + 1) * (twoNZ + 1));
                if constexpr (dim == 2) {
                    nNode = 2 * twoNX + 2 * twoNY;
                } else if constexpr (dim == 3) {
                    nNode = 2 * twoNX * twoNY + 2 * twoNY * twoNZ + 2 * twoNX * twoNZ + 2;
                }
                boundaryNode.reserve(nNode);
                for (int iZ = 0; iZ <= twoNZ; ++iZ) {
                    for (int iY = 0; iY <= twoNY; ++iY) {
                        for (int iX = 0; iX <= twoNX; ++iX) {
                            if ((iX == 0) || (iX == twoNX) || (iY == 0) || (iY == twoNY) || (iZ == 0) ||
                                (iZ == twoNZ)) {
                                boundaryNode.push_back(int(vList.size()));
                            }
                            coord[0] = 0.5 * iX * hX;
                            coord[1] = 0.5 * iY * hY;
                            coord[2] = 0.5 * iZ * hZ;
                            vList.push_back(coord);
                        }
                    }
                }
                break;
            }
            default: {
                vList.reserve((nX + 1) * (nY + 1) * (nZ + 1));
                if constexpr (dim == 2) {
                    nNode = 2 * nX + 2 * nY;
                } else if constexpr (dim == 3) {
                    nNode = 2 * nX * nY + 2 * nX * nZ + 2 * nY * nZ + 2;
                }
                boundaryNode.reserve(nNode);
                for (int iZ = 0; iZ <= nZ; ++iZ) {
                    for (int iY = 0; iY <= nY; ++iY) {
                        for (int iX = 0; iX <= nX; ++iX) {
                            if ((iX == 0) || (iX == nX) || (iY == 0) || (iY == nY) || (iZ == 0) || (iZ == nZ)) {
                                boundaryNode.push_back(vList.size());
                            }
                            coord[0] = params.lowerCorner[0] + iX * hX;
                            coord[1] = params.lowerCorner[1] + iY * hY;
                            coord[2] = params.lowerCorner[2] + iZ * hZ;
                            vList.push_back(coord);
                        }
                    }
                }
                break;
            }
        }
        return Mesh{dim, vList, std::move(cType), std::move(cList), std::move(boundaryNode)};
    }

    Mesh GenerateMesh(DomainParams const &params, std::vector<double> corners) {

        switch (params.omega) {
            case DomainType::InputFile: {
                /// To Be Implemented
            }
            case DomainType::Bar: {
                return GenerateUniformTensor<1>(params);
            }
            case DomainType::Trapeze: {
                /// To Be Implemented
            }
            case DomainType::Rectangle: {
                return GenerateUniformTensor<2>(params);
            }
            case DomainType::Brick: {
                return GenerateUniformTensor<3>(params);
            }
        }

    }

}
