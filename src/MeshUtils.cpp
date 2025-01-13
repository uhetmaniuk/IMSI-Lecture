#include "MeshUtils.h"

#include <iostream>

#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>
#include <KokkosKernels_Handle.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <KokkosSparse_spgemm.hpp>
#include <KokkosSparse_Utils.hpp>

namespace IMSI {

    template<unsigned int dim>
    Mesh GenerateUniformTensor(DomainParams const &params) {

        static_assert((dim > 0) && (dim <= 3));

        int const nX = params.numElePerDir[0];
        int const nY = (dim < 2) ? 0 : params.numElePerDir[1];
        int const nZ = (dim < 3) ? 0 : params.numElePerDir[2];
        int const numCells = nX * nY * nZ;
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

        std::vector<ElementType> cType;
        cType.resize(numCells, params.cellType);

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
                                boundaryNode.push_back(vList.size());
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

        return Mesh{std::move(vList), std::move(cType), std::move(cList), std::move(boundaryNode)};
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

    template<typename Device, typename Idx = int>
    Kokkos::StaticCrsGraph<Idx, Device>
    CombineGraphs(const Kokkos::StaticCrsGraph<Idx, Device> &aTob, const Kokkos::StaticCrsGraph<Idx, Device> &bToc) {
        //
        typedef Kokkos::StaticCrsGraph<Idx, Device> OutputGraph;
        using row_map_type = typename OutputGraph::row_map_type::non_const_type;
        using entries_type = typename OutputGraph::entries_type::non_const_type;
        //
        //row_map_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "non_const_lnow_row"), aTob.numRows() + 1);
        //--- row_mapC is initialized to 0 by default
        row_map_type row_mapC("RowMapC", aTob.numRows() + 1);
        Kokkos::parallel_for(aTob.numRows(), KOKKOS_LAMBDA(const int i) {
            row_mapC(i + 1) += 1;
            for (size_t j = aTob.row_map[i]; j < aTob.row_map[i + 1]; ++j) {
                auto const bIdx = aTob.entries[j];
                row_mapC(i + 1) += bToc.row_map[bIdx + 1] - bToc.row_map[bIdx] - 1;
            }
        });
        // create_mirror_view will only create a new view if the original one is not in HostSpace.
        auto h_row = Kokkos::create_mirror_view(row_mapC);
        h_row(0) = 0;
        for (int i = 0; i < aTob.numRows(); ++i) {
            h_row(i + 1) += h_row(i);
        }
        Kokkos::deep_copy(row_mapC, h_row);
        //
        row_map_type tmp_row("TmpRow", aTob.numRows() + 1);
        entries_type tmp_entries("TmpEntries", h_row(aTob.numRows()));
        //
        Kokkos::parallel_for(aTob.numRows(), KOKKOS_LAMBDA(const int ia) {
            for (size_t j = aTob.row_map[ia]; j < aTob.row_map[ia + 1]; ++j) {
                auto const bIdx = aTob.entries[j];
                for (size_t k = bToc.row_map[bIdx]; k < bToc.row_map[bIdx + 1]; ++k) {
                    auto cIdx = bToc.entries[k];
                    bool isStored = false;
                    for (size_t l = 0; l < tmp_row(ia + 1); ++l) {
                        if (tmp_entries(row_mapC(ia) + l) == cIdx) {
                            isStored = true;
                            break;
                        }
                    }
                    if (!isStored) {
                        tmp_entries(row_mapC(ia) + tmp_row(ia + 1)) = cIdx;
                        tmp_row(ia + 1) += 1;
                    }
                }
            }
        });
        //
        h_row = Kokkos::create_mirror_view(tmp_row);
        h_row(0) = 0;
        for (int i = 0; i < aTob.numRows(); ++i) {
            h_row(i + 1) += h_row(i);
        }
        Kokkos::deep_copy(tmp_row, h_row);
        //--- Compress the temporary entries array into `entriesC`
        entries_type entriesC("Entries", h_row(aTob.numRows()));
        Kokkos::parallel_for(aTob.numRows(), KOKKOS_LAMBDA(const int ia) {
            auto const len = tmp_row(ia + 1) - tmp_row(ia);
            for (size_t j = 0; j < len; ++j) {
                entriesC(tmp_row(ia) + j) = tmp_entries(row_mapC(ia) + j);
            }
        });
        for (size_t ia = 0; ia < aTob.numRows(); ++ia) {
            Kokkos::sort(entriesC, h_row(ia), h_row(ia + 1));
        }
        Kokkos::deep_copy(row_mapC, tmp_row);
        //
        return {entriesC, row_mapC};
    }

    void GetMatrixSparsity(const Mesh &grid) {
        typedef Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space> device_type;
        typedef Kokkos::StaticCrsGraph<int, device_type> StaticCrsGraphType;

        // Get the element-to-node in Kokkos format
        auto e2n = Kokkos::create_staticcrsgraph<StaticCrsGraphType>("CellToNode", grid.CellToNode());

        // Make the node-to-element connectivity in Kokkos format
        Kokkos::View<size_t *, Kokkos::DefaultHostExecutionSpace> rowPtr("TransposeGraphRow",
                                                                         grid.NumberVertices() + 1);
        Kokkos::View<int *, Kokkos::DefaultHostExecutionSpace> entries("TransposeGraphEntries",
                                                                       e2n.row_map[e2n.numRows()]);
        KokkosSparse::Impl::transpose_graph<StaticCrsGraphType::row_map_type,
                StaticCrsGraphType::entries_type,
                StaticCrsGraphType::row_map_type::non_const_type,
                StaticCrsGraphType::entries_type,
                StaticCrsGraphType::row_map_type::non_const_type, Kokkos::DefaultHostExecutionSpace>(grid.NumberCells(),
                                                                                                     grid.NumberVertices(),
                                                                                                     e2n.row_map,
                                                                                                     e2n.entries,
                                                                                                     rowPtr,
                                                                                                     entries);

        StaticCrsGraphType n2e(entries, rowPtr);

        // Make the node-to-node connectivity
        auto n2n = CombineGraphs<device_type, int>(n2e, e2n);
        for (int ii = 0; ii < n2n.numRows(); ++ii) {
            std::cout << ii << " = ";
            for (int j = n2n.row_map[ii]; j < n2n.row_map[ii + 1]; ++j) {
                std::cout << n2n.entries[j] << " ";
            }
            std::cout << "\n";
        }

        // Make the cell-to-cell connectivity
        auto e2e = CombineGraphs<device_type, int>(e2n, n2e);
        for (int ii = 0; ii < e2e.numRows(); ++ii) {
            std::cout << ii << " # ";
            for (int j = e2e.row_map[ii]; j < e2e.row_map[ii + 1]; ++j) {
                std::cout << e2e.entries[j] << " ";
            }
            std::cout << "\n";
        }

    }

}
