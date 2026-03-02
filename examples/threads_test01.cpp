// threads_test01.cpp
// Kokkos C++17 translation of threads_test01.jl
//
// Implements Q1/bilinear FEM stiffness assembly on a uniform nx×ny mesh
// using three strategies that mirror the Julia original:
//   1. assemble_serial   ← Julia's assemble_scalar!
//   2. assemble_threads  ← Julia's assemble_threads! (@threads)
//   3. assemble_spawn    ← Julia's assemble_spawn!  (@spawn / @sync)
//
// The key parallelism-safety mechanism (4-color checkerboard) is unchanged:
// within a color group no two elements share a node, so parallel writes to
// the sparse matrix value array are race-free.
//
// Build (OpenMP backend, recommended):
//   cmake -DCMAKE_BUILD_TYPE=Release \
//         -DKokkos_ENABLE_OPENMP=ON  \
//         -DKokkos_ARCH_NATIVE=ON    \
//         -S . -B build && cmake --build build -j
//
// Run with multiple threads:
//   OMP_NUM_THREADS=8 ./threads_test01

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

// =============================================================================
// Tiny fixed-size helpers (replace StaticArrays / SMatrix / SVector)
// =============================================================================

// 4-vector of doubles
struct Vec4 {
    double v[4]{};
    KOKKOS_INLINE_FUNCTION double  operator[](int i) const { return v[i]; }
    KOKKOS_INLINE_FUNCTION double& operator[](int i)       { return v[i]; }
};

// 2×4 matrix (column-major like Julia: col j, row i → data[j*2+i])
struct Mat24 {
    double data[8]{};
    KOKKOS_INLINE_FUNCTION double  operator()(int i, int j) const { return data[j*2+i]; }
    KOKKOS_INLINE_FUNCTION double& operator()(int i, int j)       { return data[j*2+i]; }
};

// 4×4 symmetric element stiffness matrix (row-major storage)
struct Mat44 {
    double data[16]{};
    KOKKOS_INLINE_FUNCTION double  operator()(int i, int j) const { return data[i*4+j]; }
    KOKKOS_INLINE_FUNCTION double& operator()(int i, int j)       { return data[i*4+j]; }
    KOKKOS_INLINE_FUNCTION void zero() { for (int k = 0; k < 16; ++k) data[k] = 0.0; }
    KOKKOS_INLINE_FUNCTION void add(const Mat44& o) {
        for (int k = 0; k < 16; ++k) data[k] += o.data[k];
    }
};

// =============================================================================
// Mesh data (host-side plain structs; Views for device-side kernel access)
// =============================================================================

struct Mesh {
    int nnodes, nelems, nx, ny;

    // host storage
    std::vector<double> coords;  // [2*nnodes]: (x,y) of node n → coords[2n], coords[2n+1]
    std::vector<int>    conn;    // [4*nelems]: 4 node indices (0-based) of element e
    std::vector<int>    colors;  // [nelems]: color id 0..3
    std::vector<std::vector<int>> color_groups; // 4 groups of element indices

    // Kokkos Views (1-D, zero-based) for device-accessible kernels
    Kokkos::View<double*> d_coords;   // length 2*nnodes
    Kokkos::View<int*>    d_conn;     // length 4*nelems
};

// =============================================================================
// make_mesh — mirrors make_mesh(nx, ny) in Julia (0-based indices internally)
// =============================================================================
Mesh make_mesh(int nx, int ny)
{
    Mesh m;
    m.nx     = nx;
    m.ny     = ny;
    m.nnodes = (nx + 1) * (ny + 1);
    m.nelems = nx * ny;

    m.coords.resize(2 * m.nnodes);
    m.conn.resize(4 * m.nelems);
    m.colors.resize(m.nelems);

    // Node coordinates
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            int n = j * (nx + 1) + i;          // 0-based node index
            m.coords[2*n]   = (double)i / nx;
            m.coords[2*n+1] = (double)j / ny;
        }
    }

    // Connectivity + 4-color checkerboard coloring
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int sw = j * (nx + 1) + i;          // SW corner (0-based)
            int e  = j * nx + i;                // element index (0-based)
            m.conn[4*e + 0] = sw;               // SW
            m.conn[4*e + 1] = sw + 1;           // SE
            m.conn[4*e + 2] = sw + (nx+1) + 1; // NE
            m.conn[4*e + 3] = sw + (nx+1);      // NW
            m.colors[e]     = (i & 1) + 2*(j & 1); // 0..3 (Julia used 1..4)
        }
    }

    // Color groups
    m.color_groups.resize(4);
    for (int e = 0; e < m.nelems; ++e)
        m.color_groups[m.colors[e]].push_back(e);

    // Copy to Kokkos Views
    m.d_coords = Kokkos::View<double*>("coords", 2 * m.nnodes);
    m.d_conn   = Kokkos::View<int*>  ("conn",   4 * m.nelems);
    auto hc = Kokkos::create_mirror_view(m.d_coords);
    auto hk = Kokkos::create_mirror_view(m.d_conn);
    for (int i = 0; i < 2 * m.nnodes; ++i) hc(i) = m.coords[i];
    for (int i = 0; i < 4 * m.nelems; ++i) hk(i) = m.conn[i];
    Kokkos::deep_copy(m.d_coords, hc);
    Kokkos::deep_copy(m.d_conn,   hk);

    return m;
}

// =============================================================================
// element_Ke — mirrors element_Ke(xy) in Julia
//   Computes the 4×4 Laplace stiffness matrix for one Q1 element using
//   2×2 Gauss quadrature.  KOKKOS_INLINE_FUNCTION → callable from device.
// =============================================================================
KOKKOS_INLINE_FUNCTION
Mat44 element_Ke(const Mat24& xy)
{
    const double g = 1.0 / Kokkos::sqrt(3.0);
    Mat44 Ke;
    Ke.zero();

    const double gp[2] = {-g, g};

    for (int qi = 0; qi < 2; ++qi) {
        for (int qj = 0; qj < 2; ++qj) {
            const double xi  = gp[qi];
            const double eta = gp[qj];

            // Shape function derivatives in reference space (local node 0..3 = SW,SE,NE,NW)
            Vec4 dNdxi, dNdeta;
            dNdxi[0]  = -(1.0 - eta) * 0.25;
            dNdxi[1]  =  (1.0 - eta) * 0.25;
            dNdxi[2]  =  (1.0 + eta) * 0.25;
            dNdxi[3]  = -(1.0 + eta) * 0.25;

            dNdeta[0] = -(1.0 - xi) * 0.25;
            dNdeta[1] = -(1.0 + xi) * 0.25;
            dNdeta[2] =  (1.0 + xi) * 0.25;
            dNdeta[3] =  (1.0 - xi) * 0.25;

            // Jacobian J = xy[2x4] * [dNdxi | dNdeta][4x2]  →  J[2x2]
            double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
            for (int k = 0; k < 4; ++k) {
                J11 += xy(0,k) * dNdxi[k];
                J12 += xy(0,k) * dNdeta[k];
                J21 += xy(1,k) * dNdxi[k];
                J22 += xy(1,k) * dNdeta[k];
            }
            const double detJ     = J11*J22 - J12*J21;
            const double inv_detJ = 1.0 / detJ;

            // Physical gradients: dN/dx, dN/dy via J^{-1}
            Vec4 dNdx, dNdy;
            for (int k = 0; k < 4; ++k) {
                dNdx[k] = ( J22 * dNdxi[k] - J21 * dNdeta[k]) * inv_detJ;
                dNdy[k] = (-J12 * dNdxi[k] + J11 * dNdeta[k]) * inv_detJ;
            }

            // Ke += detJ * (dNdx * dNdx' + dNdy * dNdy')
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    Ke(a,b) += detJ * (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]);
        }
    }
    return Ke;
}

// Convenience overload: gather xy from global arrays and call above
// Mirror of the @inline wrapper in Julia.
KOKKOS_INLINE_FUNCTION
Mat44 element_Ke(const Kokkos::View<double*>& coords,
                 const Kokkos::View<int*>&    conn,
                 int e)
{
    Mat24 xy;
    for (int k = 0; k < 4; ++k) {
        int n       = conn[4*e + k];
        xy(0, k)    = coords[2*n];
        xy(1, k)    = coords[2*n + 1];
    }
    return element_Ke(xy);
}

// =============================================================================
// Sparse matrix in CSR format
//   row_ptr[0..nnodes]  — start of each row in col_idx / values
//   col_idx[nnz]        — column indices
//   values[nnz]         — nonzero values (zero-initialised for assembly)
//
// elem_ptrs[16 * nelems]: for element e and local pair (a, b) (a=row, b=col,
//   both 0..3), elem_ptrs[(b*4 + a) + 16*e] = position in values[] where
//   values[conn[a,e], conn[b,e]] lives.  Exactly mirrors the Julia CSC version
//   but with row-major CSR storage instead.
// =============================================================================
struct SparseCSR {
    int nnodes, nnz_count;
    std::vector<int>    row_ptr;
    std::vector<int>    col_idx;
    std::vector<double> values;

    Kokkos::View<int*>    d_row_ptr;
    Kokkos::View<int*>    d_col_idx;
    Kokkos::View<double*> d_values;

    void sync_to_device() {
        auto hr = Kokkos::create_mirror_view(d_row_ptr);
        auto hc = Kokkos::create_mirror_view(d_col_idx);
        auto hv = Kokkos::create_mirror_view(d_values);
        for (int i = 0; i <= nnodes; ++i) hr(i) = row_ptr[i];
        for (int i = 0; i < nnz_count; ++i) { hc(i) = col_idx[i]; hv(i) = values[i]; }
        Kokkos::deep_copy(d_row_ptr, hr);
        Kokkos::deep_copy(d_col_idx, hc);
        Kokkos::deep_copy(d_values,  hv);
    }

    void sync_from_device() {
        auto hv = Kokkos::create_mirror_view(d_values);
        Kokkos::deep_copy(hv, d_values);
        for (int i = 0; i < nnz_count; ++i) values[i] = hv(i);
    }
};

// =============================================================================
// build_sparse_K — mirrors build_sparse_K(coords, conn, nnodes)
//
//   1. Collect all (row, col) pairs from element connectivities.
//   2. Sort + deduplicate to build CSR structure.
//   3. Build elem_ptrs via binary search (same O(log nnz/row) strategy).
// =============================================================================
std::pair<SparseCSR, std::vector<int>>
build_sparse_K(const Mesh& m)
{
    const int nnodes  = m.nnodes;
    const int nelems  = m.nelems;

    // --- Collect (row, col) pairs ---
    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(16 * nelems);
    for (int e = 0; e < nelems; ++e)
        for (int b = 0; b < 4; ++b)
            for (int a = 0; a < 4; ++a)
                pairs.push_back({m.conn[4*e+a], m.conn[4*e+b]});

    // --- Sort and deduplicate → CSR sparsity pattern ---
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    const int nnz = (int)pairs.size();

    SparseCSR K;
    K.nnodes    = nnodes;
    K.nnz_count = nnz;
    K.row_ptr.resize(nnodes + 1, 0);
    K.col_idx.resize(nnz);
    K.values.resize(nnz, 0.0);

    for (auto& [r, c] : pairs) K.row_ptr[r + 1]++;
    for (int i = 0; i < nnodes; ++i) K.row_ptr[i+1] += K.row_ptr[i];
    for (int k = 0; k < nnz; ++k) K.col_idx[k] = pairs[k].second;

    // Allocate device Views
    K.d_row_ptr = Kokkos::View<int*>   ("row_ptr", nnodes + 1);
    K.d_col_idx = Kokkos::View<int*>   ("col_idx", nnz);
    K.d_values  = Kokkos::View<double*>("values",  nnz);
    K.sync_to_device();

    // --- Build elem_ptrs: (b*4 + a) + 16*e → index into values[] ---
    //     For each (row=conn[a,e], col=conn[b,e]), do binary search within row.
    std::vector<int> elem_ptrs(16 * nelems);
    for (int e = 0; e < nelems; ++e) {
        for (int b = 0; b < 4; ++b) {
            int col = m.conn[4*e + b];
            for (int a = 0; a < 4; ++a) {
                int row  = m.conn[4*e + a];
                int rstart = K.row_ptr[row];
                int rend   = K.row_ptr[row + 1];
                // Binary search for col in K.col_idx[rstart..rend)
                auto it = std::lower_bound(K.col_idx.begin() + rstart,
                                           K.col_idx.begin() + rend,
                                           col);
                assert(it != K.col_idx.begin() + rend && *it == col);
                elem_ptrs[(b*4 + a) + 16*e] = (int)(it - K.col_idx.begin());
            }
        }
    }

    return {std::move(K), std::move(elem_ptrs)};
}

// =============================================================================
// Assembly helpers — inner scatter for a single element
// =============================================================================
KOKKOS_INLINE_FUNCTION
void scatter_element(const Mat44& Ke, int e,
                     const Kokkos::View<int*>&    elem_ptrs_v,
                     const Kokkos::View<double*>& values)
{
    for (int b = 0; b < 4; ++b)
        for (int a = 0; a < 4; ++a)
            values(elem_ptrs_v((b*4 + a) + 16*e)) += Ke(a, b);
}

// =============================================================================
// assemble_serial — mirrors assemble_scalar! in Julia
// =============================================================================
void assemble_serial(SparseCSR& K,
                     const Mesh& m,
                     const Kokkos::View<int*>& elem_ptrs_v,
                     const std::vector<std::vector<int>>& color_groups)
{
    Kokkos::deep_copy(K.d_values, 0.0);

    for (const auto& group : color_groups) {
        for (int e : group) {
            Mat44 Ke = element_Ke(m.d_coords, m.d_conn, e);
            scatter_element(Ke, e, elem_ptrs_v, K.d_values);
        }
    }
}

// =============================================================================
// assemble_threads — mirrors assemble_threads! (@threads) in Julia
//
//   Kokkos::parallel_for maps to the active execution space (e.g. OpenMP).
//   The barrier between color groups is implicit: parallel_for + fence()
//   both complete before the next group's loop begins, exactly like Julia's
//   implicit barrier at the end of @threads.
//
//   Race-freedom guarantee: same as Julia — within one color group, no two
//   elements share a node, so elem_ptrs positions are disjoint across threads.
// =============================================================================
void assemble_threads(SparseCSR& K,
                      const Mesh& m,
                      const Kokkos::View<int*>&    elem_ptrs_v,
                      const Kokkos::View<int*>&    group_elems_v, // flat array
                      const std::vector<int>&       group_offsets, // per-color start
                      const std::vector<std::vector<int>>& color_groups)
{
    Kokkos::deep_copy(K.d_values, 0.0);

    const auto d_coords    = m.d_coords;
    const auto d_conn      = m.d_conn;
    const auto d_values    = K.d_values;

    for (int c = 0; c < 4; ++c) {
        const int offset = group_offsets[c];
        const int count  = (int)color_groups[c].size();

        Kokkos::parallel_for(
            "assemble_color_" + std::to_string(c),
            Kokkos::RangePolicy<>(0, count),
            KOKKOS_LAMBDA(const int idx) {
                const int e  = group_elems_v(offset + idx);
                Mat44     Ke = element_Ke(d_coords, d_conn, e);
                scatter_element(Ke, e, elem_ptrs_v, d_values);
            }
        );
        // Implicit barrier: parallel_for returns only after all threads finish
        // (mirrors the end-of-@threads barrier in Julia).
        Kokkos::fence();
    }
}

// =============================================================================
// assemble_spawn — mirrors assemble_spawn! (@spawn / @sync) in Julia
//
//   Julia's @spawn partitions each color group into nthreads() chunks and
//   submits one task per chunk, with @sync as the group barrier.  The intent
//   is work-stealing composability rather than raw throughput.
//
//   Kokkos equivalent: TeamPolicy with one team per chunk, each team
//   processing a contiguous slice of the group.  This exposes the same
//   chunked structure while letting Kokkos choose scheduling.  The fence()
//   after each color provides the @sync barrier.
// =============================================================================
void assemble_spawn(SparseCSR& K,
                    const Mesh& m,
                    const Kokkos::View<int*>&    elem_ptrs_v,
                    const Kokkos::View<int*>&    group_elems_v,
                    const std::vector<int>&       group_offsets,
                    const std::vector<std::vector<int>>& color_groups)
{
    Kokkos::deep_copy(K.d_values, 0.0);

    using ExecSpace  = Kokkos::DefaultExecutionSpace;
    using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicy::member_type;

    const int nteams = ExecSpace::concurrency(); // analogous to nthreads()

    const auto d_coords = m.d_coords;
    const auto d_conn   = m.d_conn;
    const auto d_values = K.d_values;

    for (int c = 0; c < 4; ++c) {
        const int offset = group_offsets[c];
        const int ne     = (int)color_groups[c].size();
        // chunk size — mirrors cld(ne, nt) in Julia
        const int chunk  = (ne + nteams - 1) / nteams;

        Kokkos::parallel_for(
            "spawn_color_" + std::to_string(c),
            TeamPolicy(nteams, 1),          // nteams teams, 1 thread each (serial team)
            KOKKOS_LAMBDA(const TeamMember& team) {
                const int t      = team.league_rank();
                const int istart = t * chunk;
                const int iend   = Kokkos::min(istart + chunk, ne);
                // Each team (= task in Julia) processes its slice independently
                for (int idx = istart; idx < iend; ++idx) {
                    const int e  = group_elems_v(offset + idx);
                    Mat44     Ke = element_Ke(d_coords, d_conn, e);
                    scatter_element(Ke, e, elem_ptrs_v, d_values);
                }
            }
        );
        // @sync barrier — all teams must finish before next color
        Kokkos::fence();
    }
}

// =============================================================================
// Timing helper
// =============================================================================
template <class F>
double time_ms(F&& fn, int repeats = 5)
{
    fn(); Kokkos::fence(); // warm-up
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < repeats; ++r) { fn(); Kokkos::fence(); }
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / repeats;
}

// =============================================================================
// Verification helpers
// =============================================================================
double max_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double mx = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        mx = std::max(mx, std::abs(a[i] - b[i]));
    return mx;
}

// Row sums of assembled K (using CSR host arrays)
// Interior nodes of a Laplacian assembled without BCs should have row-sum ≈ 0
double max_interior_row_sum(const SparseCSR& K, const Mesh& m)
{
    const int nx = m.nx, ny = m.ny;
    double mx = 0.0;
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            int n  = j*(nx+1) + i;  // 0-based
            double s = 0.0;
            for (int k = K.row_ptr[n]; k < K.row_ptr[n+1]; ++k)
                s += K.values[k];
            mx = std::max(mx, std::abs(s));
        }
    }
    return mx;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
        const int nx = 512, ny = 512;

        // ------------------------------------------------------------------
        printf("Building mesh (%dx%d) …\n", nx, ny);
        Mesh m = make_mesh(nx, ny);
        printf("  Nodes: %d   Elements: %d   Colors: 4\n", m.nnodes, m.nelems);
        printf("  Execution space: %s\n", Kokkos::DefaultExecutionSpace::name());
        printf("  Concurrency: %d\n\n", (int)Kokkos::DefaultExecutionSpace::concurrency());

        // ------------------------------------------------------------------
        printf("Building sparsity pattern and elem_ptrs …\n");
        auto t_build0 = std::chrono::high_resolution_clock::now();
        auto [K_template, elem_ptrs_host] = build_sparse_K(m);
        auto t_build1 = std::chrono::high_resolution_clock::now();
        double build_ms = std::chrono::duration<double,std::milli>(t_build1-t_build0).count();
        printf("  K: %dx%d, %d non-zeros  [%.1f ms]\n\n",
               m.nnodes, m.nnodes, K_template.nnz_count, build_ms);

        // Upload elem_ptrs to device
        Kokkos::View<int*> elem_ptrs_v("elem_ptrs", (int)elem_ptrs_host.size());
        {
            auto h = Kokkos::create_mirror_view(elem_ptrs_v);
            for (int i = 0; i < (int)elem_ptrs_host.size(); ++i) h(i) = elem_ptrs_host[i];
            Kokkos::deep_copy(elem_ptrs_v, h);
        }

        // Flat device array of element indices per color + offsets
        // (avoids passing std::vector into KOKKOS_LAMBDA)
        std::vector<int> group_offsets(5, 0);
        for (int c = 0; c < 4; ++c)
            group_offsets[c+1] = group_offsets[c] + (int)m.color_groups[c].size();
        const int total_elems = group_offsets[4];

        Kokkos::View<int*> group_elems_v("group_elems", total_elems);
        {
            auto h = Kokkos::create_mirror_view(group_elems_v);
            for (int c = 0; c < 4; ++c)
                for (int i = 0; i < (int)m.color_groups[c].size(); ++i)
                    h(group_offsets[c] + i) = m.color_groups[c][i];
            Kokkos::deep_copy(group_elems_v, h);
        }

        // Three separate K copies (mirrors K1/K2/K3 in Julia)
        SparseCSR K1 = K_template;
        SparseCSR K2 = K_template;
        SparseCSR K3 = K_template;
        // Each needs its own device values view
        K1.d_values = Kokkos::View<double*>("values1", K_template.nnz_count);
        K2.d_values = Kokkos::View<double*>("values2", K_template.nnz_count);
        K3.d_values = Kokkos::View<double*>("values3", K_template.nnz_count);

        // ------------------------------------------------------------------
        // Assemble once each to verify correctness
        assemble_serial (K1, m, elem_ptrs_v, m.color_groups);
        assemble_threads(K2, m, elem_ptrs_v, group_elems_v, group_offsets, m.color_groups);
        assemble_spawn  (K3, m, elem_ptrs_v, group_elems_v, group_offsets, m.color_groups);

        K1.sync_from_device();
        K2.sync_from_device();
        K3.sync_from_device();

        printf("Max error (RangePolicy) : %.3e\n", max_diff(K1.values, K2.values));
        printf("Max error (TeamPolicy)   : %.3e\n", max_diff(K1.values, K3.values));
        printf("Max |row-sum| on interior nodes: %.3e\n\n",
               max_interior_row_sum(K1, m));

        // ------------------------------------------------------------------
        // Benchmarks (mirrors @btime in Julia)
        const int reps = 512;
        printf("Serial   : %.3f ms\n",
               time_ms([&]{ assemble_serial (K1, m, elem_ptrs_v, m.color_groups); }, reps));
        printf("RangePolicy : %.3f ms\n",
               time_ms([&]{ assemble_threads(K2, m, elem_ptrs_v, group_elems_v, group_offsets, m.color_groups); }, reps));
        printf("TeamPolicy  : %.3f ms\n",
               time_ms([&]{ assemble_spawn  (K3, m, elem_ptrs_v, group_elems_v, group_offsets, m.color_groups); }, reps));
    }
    Kokkos::finalize();
    return 0;
}

