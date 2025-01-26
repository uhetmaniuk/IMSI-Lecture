#pragma once

#include "Element.h"
#include "Mesh.h"

#include <Kokkos_StaticCrsGraph.hpp>

#include <string>
#include <vector>

namespace IMSI {

    enum class DomainType : char {
        InputFile,
        Bar,
        Rectangle, Trapeze,
        Brick
    };


    struct DomainParams {
        std::string fileName;
        double lowerCorner[3] = {0.0, 0.0, 0.0};
        double upperCorner[3] = {1.0, 1.0, 1.0};
        int numElePerDir[3] = {0, 0, 0};
        DomainType omega = DomainType::Rectangle;
        ElementType cellType = ElementType::Q1;
    };

    template<typename Device = Kokkos::DefaultHostExecutionSpace>
    struct MeshConnectivity {

        /// \brief Reference to mesh
        const Mesh &mesh;

        /// \brief Node-to-node connectivity
        Kokkos::StaticCrsGraph<int, Device> n2n;

        /// \brief Element-to-element connectivity
        Kokkos::StaticCrsGraph<int, Device> e2e;

        /// \brief Color to element connectivity
        Kokkos::StaticCrsGraph<int, Device> c2e;
    };

    /// \brief Function to generate the finite element mesh
    Mesh GenerateMesh(DomainParams const &params, std::vector<double> corners = {});

    template<typename Device, typename Idx = int>
    Kokkos::StaticCrsGraph <Idx, Device>
    CombineGraphs(const Kokkos::StaticCrsGraph <Idx, Device> &aTob, const Kokkos::StaticCrsGraph <Idx, Device> &bToc);

    template<typename Device, typename Idx = int>
    Kokkos::StaticCrsGraph <Idx, Device>
    TransposeGraph(const Kokkos::StaticCrsGraph <Idx, Device> &aTob, int numRowsB);

    template<typename Device, typename Idx = int>
    Kokkos::StaticCrsGraph <Idx, Device>
    ColorGraph(const Kokkos::StaticCrsGraph <Idx, Device> &e2e);

    template<typename Device, typename Idx = int>
    void SortEntries(Kokkos::StaticCrsGraph <Idx, Device> &g);

    template<typename Device = Kokkos::DefaultHostExecutionSpace>
    MeshConnectivity<Device> GetMeshConnectivity(const Mesh &grid);

}

//
// Definition of templated functions
//

#include <Kokkos_Core.hpp>
#include <KokkosGraph_Distance1Color.hpp>
#include <KokkosKernels_Handle.hpp>
#include <KokkosSparse_Utils.hpp>

namespace IMSI {

    template<typename Device, typename Idx>
    Kokkos::StaticCrsGraph <Idx, Device>
    CombineGraphs(const Kokkos::StaticCrsGraph <Idx, Device> &aTob, const Kokkos::StaticCrsGraph <Idx, Device> &bToc) {
        typedef Kokkos::StaticCrsGraph<Idx, Device> OutputGraph;
        using row_map_type = typename OutputGraph::row_map_type::non_const_type;
        using entries_type = typename OutputGraph::entries_type::non_const_type;
        //
        //--- row_mapC is initialized to 0 by default
        //row_map_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "non_const_lnow_row"), aTob.numRows() + 1);
        row_map_type row_mapC("RowMapC", aTob.numRows() + 1);
        Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, aTob.numRows()),
                             KOKKOS_LAMBDA(const int i) {
                                 row_mapC(i + 1) += 1;
                                 for (size_t j = aTob.row_map[i]; j < aTob.row_map[i + 1]; ++j) {
                                     auto const bIdx = aTob.entries[j];
                                     row_mapC(i + 1) += bToc.row_map[bIdx + 1] - bToc.row_map[bIdx] - 1;
                                 }
                             });
        // create_mirror_view will only create a new view if the original one is not in HostSpace.
        auto h_row = Kokkos::create_mirror_view(row_mapC);
        Kokkos::deep_copy(h_row, row_mapC);
        h_row(0) = 0;
        for (int i = 0; i < aTob.numRows(); ++i) {
            h_row(i + 1) += h_row(i);
        }
        Kokkos::deep_copy(row_mapC, h_row);
        //
        row_map_type tmp_row("TmpRow", aTob.numRows() + 1);
        entries_type tmp_entries("TmpEntries", h_row(aTob.numRows()));
        //
        Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, aTob.numRows()),
                             KOKKOS_LAMBDA(const int ia) {
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
        h_row = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), tmp_row);
        h_row(0) = 0;
        for (int i = 0; i < aTob.numRows(); ++i) {
            h_row(i + 1) += h_row(i);
        }
        Kokkos::deep_copy(tmp_row, h_row);
        //--- Compress the temporary entries array into `entriesC`
        entries_type entriesC("Entries", h_row(aTob.numRows()));
        Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, aTob.numRows()),
                             KOKKOS_LAMBDA(const int ia) {
                                 auto const len = tmp_row(ia + 1) - tmp_row(ia);
                                 for (size_t j = 0; j < len; ++j) {
                                     entriesC(tmp_row(ia) + j) = tmp_entries(row_mapC(ia) + j);
                                 }
                             });
        Kokkos::deep_copy(row_mapC, tmp_row);
        //
        return {entriesC, row_mapC};
    }

    template<typename Device, typename Idx>
    Kokkos::StaticCrsGraph <Idx, Device>
    TransposeGraph(const Kokkos::StaticCrsGraph <Idx, Device> &aTob, int numRowsB) {
        typedef Kokkos::StaticCrsGraph<Idx, Device> OutputGraph;
        using row_map_type = typename OutputGraph::row_map_type::non_const_type;
        using entries_type = typename OutputGraph::entries_type::non_const_type;
        //
        row_map_type rowPtr("TransposeGraphRow", numRowsB + 1);
        entries_type entries("TransposeGraphEntries", aTob.entries.size());
        KokkosSparse::Impl::transpose_graph<typename OutputGraph::row_map_type,
                typename OutputGraph::entries_type, row_map_type, entries_type, row_map_type, typename Device::execution_space>(
                aTob.numRows(),
                numRowsB,
                aTob.row_map,
                aTob.entries,
                rowPtr,
                entries);
        return {entries, rowPtr};
    }

    template<typename Device, typename Idx>
    void SortEntries(Kokkos::StaticCrsGraph <Idx, Device> &g) {
        auto h_row = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), g.row_map);
        auto h_entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), g.entries);
        // Sort on the host
        // TO DO: Check why Kokkos::sort is expensive
        Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, g.numRows()),
                             KOKKOS_LAMBDA(int i) {
                                 std::sort(&h_entries(h_row(i)), &h_entries(h_row(i + 1)));
                             });
        Kokkos::deep_copy(g.entries, h_entries);
    }

    template<typename Device, typename Idx>
    Kokkos::StaticCrsGraph <Idx, Device>
    ColorGraph(const Kokkos::StaticCrsGraph <Idx, Device> &e2e) {
        typedef Kokkos::StaticCrsGraph<Idx, Device> OutputGraph;
        using row_map_type = typename OutputGraph::row_map_type::non_const_type;
        using entries_type = typename OutputGraph::entries_type::non_const_type;
        //
        using DeviceSpace = typename Device::memory_space;
        KokkosKernels::Experimental::KokkosKernelsHandle<typename OutputGraph::size_type, Idx, KokkosKernels::default_scalar, Device, DeviceSpace, DeviceSpace> handle;
        // Use the default algorithm (chosen based on ExecSpace)
        handle.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
        // Run coloring
        KokkosGraph::Experimental::graph_color(&handle, e2e.numRows(), e2e.numRows(), e2e.row_map, e2e.entries);
        // Get the colors array, and the number of colors used from the handle.
        auto colors = handle.get_graph_coloring_handle()->get_vertex_colors();
        auto numColors = handle.get_graph_coloring_handle()->get_num_colors();
        // colors maps to 1-based color index
        Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, colors.size()),
                             KOKKOS_LAMBDA(const int ie) {
                                 colors(ie) -= 1;
                             });
        // Clean up
        handle.destroy_graph_coloring_handle();
        //
        row_map_type e2cPtr("EleToColorPtr", e2e.numRows() + 1);
        auto h_ptr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e2cPtr);
        for (int ii = 0; ii + 1 < h_ptr.size(); ++ii) {
            h_ptr(ii + 1) = h_ptr(ii) + 1;
        }
        Kokkos::deep_copy(e2cPtr, h_ptr);
        // Get the element-to-color in Kokkos format
        Kokkos::StaticCrsGraph<Idx, Device> e2c(colors, e2cPtr);
        //
        return TransposeGraph(e2c, numColors);
    }

    template<typename Device>
    MeshConnectivity<Device> GetMeshConnectivity(const Mesh &grid) {
        typedef Kokkos::StaticCrsGraph<int, Device> StaticCrsGraphType;

        // Get the element-to-node in Kokkos format
        auto e2n = Kokkos::create_staticcrsgraph<StaticCrsGraphType>("CellToNode", grid.CellToNode());

        // Make the node-to-element connectivity in Kokkos format
        auto start = std::chrono::high_resolution_clock::now();
        StaticCrsGraphType n2e = TransposeGraph<Device, int>(e2n, int(grid.NumberVertices()));
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = end - start;
        std::cout << " %%% n2e = " << dt.count() << "\n";

        // Make the node-to-node connectivity
        start = std::chrono::high_resolution_clock::now();
        auto n2n = CombineGraphs<Device, int>(n2e, e2n);
        SortEntries(n2n);
        end = std::chrono::high_resolution_clock::now();
        dt = end - start;
        std::cout << " %%% n2n = " << dt.count() << "\n";

        // Make the cell-to-cell connectivity
        start = std::chrono::high_resolution_clock::now();
        auto e2e = CombineGraphs<Device, int>(e2n, n2e);
        end = std::chrono::high_resolution_clock::now();
        dt = end - start;
        std::cout << " %%% e2e = " << dt.count() << "\n";

        // Coloring
        start = std::chrono::high_resolution_clock::now();
        auto c2e = ColorGraph(e2e);
        SortEntries(c2e);
        end = std::chrono::high_resolution_clock::now();
        dt = end - start;
        std::cout << " %%% c2e = " << dt.count() << "\n";

        return {grid, n2n, e2e, c2e};
    }

}
