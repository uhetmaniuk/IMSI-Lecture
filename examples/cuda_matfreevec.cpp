// cuda_matfreevec.cpp
// Kokkos C++17 translation of cuda_matfreevec.jl
//
// Implements matrix-free Q1/bilinear FEM Laplace SpMV (y = K·x) in four
// flavours that mirror the Julia original:
//
//   1. matvec_cpu_serial   ← matvec_cpu!          (host, single-thread)
//   2. matvec_gpu_atomic   ← matvec_gpu_atomic!    (one thread/elem, atomicAdd)
//   3. matvec_gpu_color    ← matvec_gpu_color!     (4-color, race-free writes)
//   4. assemble_K_cpu / SpMV                       (reference)
//
// Design notes
// ─────────────
// • Kokkos::atomic_add replaces CUDA.@atomic
// • parallel_for + RangePolicy replaces @cuda threads=256 blocks=N
// • Kokkos::fence() replaces CUDA.@sync
// • conn stored as int32_t (halves bandwidth on GPU, same as Julia Int32)
// • CSR reference matrix built host-side from sorted (row,col) pairs
//
// Build (CUDA backend):
//   cmake -DCMAKE_BUILD_TYPE=Release  \
//         -DKokkos_ENABLE_CUDA=ON     \
//         -DKokkos_ARCH_AMPERE80=ON   \
//         -S . -B build && cmake --build build -j
//
// Build (OpenMP backend, CPU-only benchmarks):
//   cmake -DCMAKE_BUILD_TYPE=Release  \
//         -DKokkos_ENABLE_OPENMP=ON   \
//         -DKokkos_ARCH_NATIVE=ON     \
//         -S . -B build && cmake --build build -j

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Convenience aliases
// ─────────────────────────────────────────────────────────────────────────────
using ExecSpace  = Kokkos::DefaultExecutionSpace;
using MemSpace   = typename ExecSpace::memory_space;
using HostExec   = Kokkos::DefaultHostExecutionSpace;
using HostMem    = typename HostExec::memory_space;

// Device Views
using DView1d   = Kokkos::View<double*,   MemSpace>;
using DView2d   = Kokkos::View<double**,  Kokkos::LayoutLeft, MemSpace>; // column-major
using DViewI32  = Kokkos::View<int32_t**, Kokkos::LayoutLeft, MemSpace>; // conn (col-major)
using DViewIdx  = Kokkos::View<int32_t*,  MemSpace>;                     // color group indices

// Host mirrors
using HView1d   = typename DView1d::HostMirror;
using HView2d   = typename DView2d::HostMirror;
using HViewI32  = typename DViewI32::HostMirror;

// ─────────────────────────────────────────────────────────────────────────────
// Mesh
// ─────────────────────────────────────────────────────────────────────────────
struct Mesh {
    int nx, ny, nnodes, nelems;

    // Host arrays
    std::vector<double>  coords_h;   // [2 * nnodes], column-major: (x,y) of node n
    std::vector<int32_t> conn_h;     // [4 * nelems], column-major: 4 nodes of elem e
    std::vector<int>     colors_h;   // [nelems]
    std::vector<std::vector<int32_t>> color_groups_h;  // 4 groups

    // Device Views
    DView2d  d_coords;   // (2, nnodes) col-major
    DViewI32 d_conn;     // (4, nelems) col-major, int32_t
};

// make_mesh — mirrors make_mesh(nx, ny) in Julia
Mesh make_mesh(int nx, int ny)
{
    Mesh m;
    m.nx     = nx;  m.ny     = ny;
    m.nnodes = (nx + 1) * (ny + 1);
    m.nelems = nx * ny;

    m.coords_h.resize(2 * m.nnodes);
    m.conn_h.resize(4 * m.nelems);
    m.colors_h.resize(m.nelems);
    m.color_groups_h.resize(4);

    // Node coordinates  (col-major: coords[2*n], coords[2*n+1])
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            int n = j * (nx + 1) + i;
            m.coords_h[2*n]     = (double)i / nx;
            m.coords_h[2*n + 1] = (double)j / ny;
        }

    // Connectivity + 4-color checkerboard
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int sw = j * (nx + 1) + i;
            int e  = j * nx + i;
            m.conn_h[4*e + 0] = (int32_t)sw;
            m.conn_h[4*e + 1] = (int32_t)(sw + 1);
            m.conn_h[4*e + 2] = (int32_t)(sw + (nx+1) + 1);
            m.conn_h[4*e + 3] = (int32_t)(sw + (nx+1));
            int c = (i & 1) + 2 * (j & 1);   // 0..3
            m.colors_h[e] = c;
            m.color_groups_h[c].push_back((int32_t)e);
        }

    // ── Upload to device ──────────────────────────────────────────────────────
    m.d_coords = DView2d("coords", 2, m.nnodes);
    m.d_conn   = DViewI32("conn",  4, m.nelems);

    auto hc = Kokkos::create_mirror_view(m.d_coords);
    auto hk = Kokkos::create_mirror_view(m.d_conn);
    for (int n = 0; n < m.nnodes; ++n) {
        hc(0, n) = m.coords_h[2*n];
        hc(1, n) = m.coords_h[2*n + 1];
    }
    for (int e = 0; e < m.nelems; ++e)
        for (int k = 0; k < 4; ++k)
            hk(k, e) = m.conn_h[4*e + k];
    Kokkos::deep_copy(m.d_coords, hc);
    Kokkos::deep_copy(m.d_conn,   hk);

    return m;
}

// ─────────────────────────────────────────────────────────────────────────────
// Gauss-point kernel — shared by all four variants
//
// Given element coordinates and local u values, accumulates the
// matrix-free contributions B'*(B*u_loc) at all 4 Gauss points.
// KOKKOS_INLINE_FUNCTION makes it callable from both device and host.
//
// Mirrors the inner for η in (-g,g), ξ in (-g,g) loop in Julia.
// ─────────────────────────────────────────────────────────────────────────────
KOKKOS_INLINE_FUNCTION
void gauss_matvec(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4,
    double u1, double u2, double u3, double u4,
    double& v1, double& v2, double& v3, double& v4)
{
    const double g = 1.0 / Kokkos::sqrt(3.0);
    const double gpts[2] = {-g, g};

    v1 = v2 = v3 = v4 = 0.0;

    for (int qi = 0; qi < 2; ++qi) {
        for (int qj = 0; qj < 2; ++qj) {
            const double xi  = gpts[qi];
            const double eta = gpts[qj];

            // Shape function derivatives in reference space
            const double dN1dxi = -(1.0 - eta) * 0.25;
            const double dN2dxi =  (1.0 - eta) * 0.25;
            const double dN3dxi =  (1.0 + eta) * 0.25;
            const double dN4dxi = -(1.0 + eta) * 0.25;

            const double dN1deta = -(1.0 - xi) * 0.25;
            const double dN2deta = -(1.0 + xi) * 0.25;
            const double dN3deta =  (1.0 + xi) * 0.25;
            const double dN4deta =  (1.0 - xi) * 0.25;

            // Jacobian
            const double J11 = x1*dN1dxi  + x2*dN2dxi  + x3*dN3dxi  + x4*dN4dxi;
            const double J12 = x1*dN1deta + x2*dN2deta + x3*dN3deta + x4*dN4deta;
            const double J21 = y1*dN1dxi  + y2*dN2dxi  + y3*dN3dxi  + y4*dN4dxi;
            const double J22 = y1*dN1deta + y2*dN2deta + y3*dN3deta + y4*dN4deta;

            const double detJ = J11*J22 - J12*J21;
            const double iJ   = 1.0 / detJ;

            // Physical gradients
            const double dN1dx = ( J22*dN1dxi - J21*dN1deta) * iJ;
            const double dN1dy = (-J12*dN1dxi + J11*dN1deta) * iJ;
            const double dN2dx = ( J22*dN2dxi - J21*dN2deta) * iJ;
            const double dN2dy = (-J12*dN2dxi + J11*dN2deta) * iJ;
            const double dN3dx = ( J22*dN3dxi - J21*dN3deta) * iJ;
            const double dN3dy = (-J12*dN3dxi + J11*dN3deta) * iJ;
            const double dN4dx = ( J22*dN4dxi - J21*dN4deta) * iJ;
            const double dN4dy = (-J12*dN4dxi + J11*dN4deta) * iJ;

            // B·u: project u onto gradient space at this Gauss point
            const double Bu_x = dN1dx*u1 + dN2dx*u2 + dN3dx*u3 + dN4dx*u4;
            const double Bu_y = dN1dy*u1 + dN2dy*u2 + dN3dy*u3 + dN4dy*u4;

            // B'*(B·u): accumulate local result, weighted by detJ (Gauss weight = 1)
            v1 += detJ * (dN1dx*Bu_x + dN1dy*Bu_y);
            v2 += detJ * (dN2dx*Bu_x + dN2dy*Bu_y);
            v3 += detJ * (dN3dx*Bu_x + dN3dy*Bu_y);
            v4 += detJ * (dN4dx*Bu_x + dN4dy*Bu_y);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 1.  matvec_cpu_serial — mirrors matvec_cpu! in Julia
//
//     Single-thread host baseline.  Iterates all elements sequentially,
//     computes B'*B*u on the fly, scatters with plain += (no conflict).
// ─────────────────────────────────────────────────────────────────────────────
void matvec_cpu_serial(
    std::vector<double>&       v,
    const std::vector<double>& u,
    const Mesh&                m)
{
    std::fill(v.begin(), v.end(), 0.0);

    const auto& coords = m.coords_h;
    const auto& conn   = m.conn_h;

    for (int e = 0; e < m.nelems; ++e) {
        const int n1 = conn[4*e+0], n2 = conn[4*e+1],
                  n3 = conn[4*e+2], n4 = conn[4*e+3];

        double v1, v2, v3, v4;
        gauss_matvec(
            coords[2*n1], coords[2*n1+1],
            coords[2*n2], coords[2*n2+1],
            coords[2*n3], coords[2*n3+1],
            coords[2*n4], coords[2*n4+1],
            u[n1], u[n2], u[n3], u[n4],
            v1, v2, v3, v4);

        v[n1] += v1;  v[n2] += v2;
        v[n3] += v3;  v[n4] += v4;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2.  matvec_gpu_atomic — mirrors matvec_gpu_atomic! in Julia
//
//     One thread per element, mapped via RangePolicy (replaces @cuda with
//     explicit blockDim/blockIdx arithmetic).  Scatter uses
//     Kokkos::atomic_add which compiles to hardware float64 atomics on
//     Volta+ (sm_70+), matching CUDA.@atomic.
// ─────────────────────────────────────────────────────────────────────────────
void matvec_gpu_atomic(
    const DView1d&   d_v,
    const DView1d&   d_u,
    const DView2d&   d_coords,
    const DViewI32&  d_conn,
    int              nelems)
{
    Kokkos::deep_copy(d_v, 0.0);

    Kokkos::parallel_for(
        "matvec_atomic",
        Kokkos::RangePolicy<ExecSpace>(0, nelems),
        KOKKOS_LAMBDA(const int e) {
            const int32_t n1 = d_conn(0, e), n2 = d_conn(1, e),
                          n3 = d_conn(2, e), n4 = d_conn(3, e);

            double v1, v2, v3, v4;
            gauss_matvec(
                d_coords(0, n1), d_coords(1, n1),
                d_coords(0, n2), d_coords(1, n2),
                d_coords(0, n3), d_coords(1, n3),
                d_coords(0, n4), d_coords(1, n4),
                d_u(n1), d_u(n2), d_u(n3), d_u(n4),
                v1, v2, v3, v4);

            // Atomic scatter — replaces CUDA.@atomic
            // Multiple elements share nodes → writes must be atomic.
            Kokkos::atomic_add(&d_v(n1), v1);
            Kokkos::atomic_add(&d_v(n2), v2);
            Kokkos::atomic_add(&d_v(n3), v3);
            Kokkos::atomic_add(&d_v(n4), v4);
        }
    );
    // Kokkos::fence() is called by the benchmark harness; kernel itself is async.
}

// ─────────────────────────────────────────────────────────────────────────────
// 3.  matvec_gpu_color — mirrors matvec_gpu_color! in Julia
//
//     Four sequential kernel launches, one per color.  Within a color group
//     no two elements share a node, so threads write directly to d_v with
//     plain += — no atomics needed.
//
//     Sequential launches in the same Kokkos execution space instance use the
//     same stream, so each parallel_for acts as an implicit barrier before the
//     next, exactly like Julia's implicit stream-ordering between @cuda calls.
//
//     Trade-off vs. atomic variant — unchanged from Julia comments:
//       + No atomic contention → better on older GPUs
//       - 4 kernel launches (small overhead per launch)
//       - Requires color-group index arrays on device
// ─────────────────────────────────────────────────────────────────────────────
void matvec_gpu_color(
    const DView1d&                    d_v,
    const DView1d&                    d_u,
    const DView2d&                    d_coords,
    const DViewI32&                   d_conn,
    const std::vector<DViewIdx>&      d_color_groups)
{
    Kokkos::deep_copy(d_v, 0.0);

    for (int c = 0; c < 4; ++c) {
        const DViewIdx& grp    = d_color_groups[c];
        const int       ngroup = (int)grp.extent(0);

        Kokkos::parallel_for(
            "matvec_color_" + std::to_string(c),
            Kokkos::RangePolicy<ExecSpace>(0, ngroup),
            KOKKOS_LAMBDA(const int idx) {
                const int32_t e  = grp(idx);
                const int32_t n1 = d_conn(0, e), n2 = d_conn(1, e),
                              n3 = d_conn(2, e), n4 = d_conn(3, e);

                double v1, v2, v3, v4;
                gauss_matvec(
                    d_coords(0, n1), d_coords(1, n1),
                    d_coords(0, n2), d_coords(1, n2),
                    d_coords(0, n3), d_coords(1, n3),
                    d_coords(0, n4), d_coords(1, n4),
                    d_u(n1), d_u(n2), d_u(n3), d_u(n4),
                    v1, v2, v3, v4);

                // Direct write — race-free because same-color elements
                // never share nodes. Mirrors the plain += in Julia.
                d_v(n1) += v1;  d_v(n2) += v2;
                d_v(n3) += v3;  d_v(n4) += v4;
            }
        );
        // Implicit stream barrier: the next parallel_for cannot start
        // until this one completes (same CUDA stream).
        // Explicit fence only needed when reading results on host.
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4.  Reference: assemble sparse K on the host + apply it
//
//     Mirrors assemble_K_cpu + K_ref * u_cpu in Julia.
//     CSR format (row-sorted).  Used only for correctness verification.
// ─────────────────────────────────────────────────────────────────────────────
struct CSR {
    int nnodes, nnz;
    std::vector<int>    row_ptr;
    std::vector<int>    col_idx;
    std::vector<double> values;
};

CSR assemble_K_cpu(const Mesh& m)
{
    // Collect all (row, col, val) triplets — 16 per element
    std::vector<std::tuple<int,int,double>> trips;
    trips.reserve(16 * m.nelems);

    const double g = 1.0 / std::sqrt(3.0);
    const double gpts[2] = {-g, g};

    for (int e = 0; e < m.nelems; ++e) {
        const int n[4] = {
            m.conn_h[4*e+0], m.conn_h[4*e+1],
            m.conn_h[4*e+2], m.conn_h[4*e+3]
        };
        const double x[4] = {
            m.coords_h[2*n[0]], m.coords_h[2*n[1]],
            m.coords_h[2*n[2]], m.coords_h[2*n[3]]
        };
        const double y[4] = {
            m.coords_h[2*n[0]+1], m.coords_h[2*n[1]+1],
            m.coords_h[2*n[2]+1], m.coords_h[2*n[3]+1]
        };

        double Ke[4][4] = {};
        for (int qi = 0; qi < 2; ++qi) {
            for (int qj = 0; qj < 2; ++qj) {
                const double xi = gpts[qi], eta = gpts[qj];
                const double dNdxi[4]  = {-(1-eta)*.25, (1-eta)*.25,  (1+eta)*.25, -(1+eta)*.25};
                const double dNdeta[4] = {-(1-xi) *.25,-(1+xi) *.25,  (1+xi) *.25,  (1-xi) *.25};
                double J11=0, J12=0, J21=0, J22=0;
                for (int k=0; k<4; ++k) {
                    J11 += x[k]*dNdxi[k];  J12 += x[k]*dNdeta[k];
                    J21 += y[k]*dNdxi[k];  J22 += y[k]*dNdeta[k];
                }
                const double detJ = J11*J22 - J12*J21, iJ = 1.0/detJ;
                double dNdx[4], dNdy[4];
                for (int k=0; k<4; ++k) {
                    dNdx[k] = ( J22*dNdxi[k] - J21*dNdeta[k]) * iJ;
                    dNdy[k] = (-J12*dNdxi[k] + J11*dNdeta[k]) * iJ;
                }
                for (int a=0; a<4; ++a)
                    for (int b=0; b<4; ++b)
                        Ke[a][b] += detJ * (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]);
            }
        }
        for (int a=0; a<4; ++a)
            for (int b=0; b<4; ++b)
                trips.emplace_back(n[a], n[b], Ke[a][b]);
    }

    // Sort by (row, col) so duplicates are adjacent
    std::sort(trips.begin(), trips.end());

    // ── Pass 1: merge duplicate (row,col) pairs, summing values ──────────────
    struct Entry { int row, col; double val; };
    std::vector<Entry> entries;
    entries.reserve(trips.size());
    for (auto& [r, c, v] : trips) {
        if (!entries.empty() && entries.back().row == r && entries.back().col == c)
            entries.back().val += v;
        else
            entries.push_back({r, c, v});
    }

    // ── Pass 2: build row_ptr via count array + exclusive prefix sum ──────────
    // row_ptr[r]   = index of first nonzero in row r
    // row_ptr[r+1] = one past last nonzero in row r  (standard CSR)
    CSR K;
    K.nnodes  = m.nnodes;
    K.nnz     = (int)entries.size();
    K.row_ptr.resize(m.nnodes + 1, 0);
    K.col_idx.resize(K.nnz);
    K.values .resize(K.nnz);

    // Count nonzeros per row into row_ptr[row+1]
    for (auto& e : entries) K.row_ptr[e.row + 1]++;

    // Exclusive prefix sum turns counts into start offsets
    for (int r = 0; r < m.nnodes; ++r) K.row_ptr[r + 1] += K.row_ptr[r];

    // Fill col_idx / values (entries already sorted by row then col)
    for (int i = 0; i < K.nnz; ++i) {
        K.col_idx[i] = entries[i].col;
        K.values[i]  = entries[i].val;
    }
    return K;
}

// y = K * x  (host CSR SpMV)
void spmv(std::vector<double>&       y,
          const CSR&                  K,
          const std::vector<double>& x)
{
    for (int r = 0; r < K.nnodes; ++r) {
        double s = 0.0;
        for (int k = K.row_ptr[r]; k < K.row_ptr[r+1]; ++k)
            s += K.values[k] * x[K.col_idx[k]];
        y[r] = s;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Timing helper
// ─────────────────────────────────────────────────────────────────────────────
template <class F>
double time_ms(F&& fn, int reps = 10)
{
    fn(); Kokkos::fence(); // warm-up
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < reps; ++r) { fn(); Kokkos::fence(); }
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / reps;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double mx = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        mx = std::max(mx, std::abs(a[i] - b[i]));
    return mx;
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    // Parse --nx= and --ny= before Kokkos::initialize so Kokkos only sees
    // its own flags.  Unknown flags are left in argv for Kokkos.
    int nx = 1024, ny = 1024;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if      (arg.rfind("--nx=", 0) == 0) nx = std::stoi(arg.substr(5));
        else if (arg.rfind("--ny=", 0) == 0) ny = std::stoi(arg.substr(5));
    }

    Kokkos::initialize(argc, argv);
    {

        // ── Mesh ─────────────────────────────────────────────────────────────
        printf("Building mesh (%d×%d)…\n", nx, ny);
        Mesh m = make_mesh(nx, ny);
        printf("  Nodes: %d   Elements: %d\n", m.nnodes, m.nelems);
        printf("  Execution space: %s\n\n", ExecSpace::name());

        // ── Random test vector ────────────────────────────────────────────────
        std::mt19937_64 rng(42);
        std::normal_distribution<double> dist(0.0, 1.0);
        std::vector<double> u_h(m.nnodes), v_ref(m.nnodes), v_cpu(m.nnodes);
        for (auto& x : u_h) x = dist(rng);

        // ── Reference: assembled sparse K (host only) ─────────────────────────
        printf("Assembling reference sparse K on CPU…\n");
        auto t0 = std::chrono::high_resolution_clock::now();
        CSR K_ref = assemble_K_cpu(m);
        auto t1 = std::chrono::high_resolution_clock::now();
        printf("  nnz = %d  [%.1f ms]\n\n",
               K_ref.nnz,
               std::chrono::duration<double,std::milli>(t1-t0).count());
        spmv(v_ref, K_ref, u_h);

        // ── CPU matrix-free ───────────────────────────────────────────────────
        matvec_cpu_serial(v_cpu, u_h, m);
        printf("Max error (CPU matrix-free)  : %.3e\n",
               max_abs_diff(v_cpu, v_ref));

        // ── Upload to device ──────────────────────────────────────────────────
        DView1d d_u("u", m.nnodes);
        DView1d d_v("v", m.nnodes);
        {
            auto hu = Kokkos::create_mirror_view(d_u);
            for (int i = 0; i < m.nnodes; ++i) hu(i) = u_h[i];
            Kokkos::deep_copy(d_u, hu);
        }

        // Color-group index arrays on device (mirrors CuArray(Int32.(g)))
        std::vector<DViewIdx> d_color_groups(4);
        for (int c = 0; c < 4; ++c) {
            const auto& g = m.color_groups_h[c];
            d_color_groups[c] = DViewIdx("color_grp_" + std::to_string(c), g.size());
            auto hg = Kokkos::create_mirror_view(d_color_groups[c]);
            for (int i = 0; i < (int)g.size(); ++i) hg(i) = g[i];
            Kokkos::deep_copy(d_color_groups[c], hg);
        }

        // ── GPU atomic variant ────────────────────────────────────────────────
        matvec_gpu_atomic(d_v, d_u, m.d_coords, m.d_conn, m.nelems);
        Kokkos::fence();
        {
            auto hv = Kokkos::create_mirror_view(d_v);
            Kokkos::deep_copy(hv, d_v);
            std::vector<double> v_gpu_atomic(m.nnodes);
            for (int i = 0; i < m.nnodes; ++i) v_gpu_atomic[i] = hv(i);
            printf("Max error (GPU atomic)       : %.3e\n",
                   max_abs_diff(v_gpu_atomic, v_ref));
        }

        // ── GPU coloring variant ──────────────────────────────────────────────
        matvec_gpu_color(d_v, d_u, m.d_coords, m.d_conn, d_color_groups);
        Kokkos::fence();
        {
            auto hv = Kokkos::create_mirror_view(d_v);
            Kokkos::deep_copy(hv, d_v);
            std::vector<double> v_gpu_color(m.nnodes);
            for (int i = 0; i < m.nnodes; ++i) v_gpu_color[i] = hv(i);
            printf("Max error (GPU coloring)     : %.3e\n",
                   max_abs_diff(v_gpu_color, v_ref));
        }
        printf("\n");

        // ── Benchmarks (mirrors @btime + CUDA.@sync in Julia) ─────────────────
        // Kokkos::fence() inside time_ms() plays the role of CUDA.@sync:
        // it blocks the host until all device work in the default stream is done.
        const int reps = 10;

        printf("CPU sparse K*u (ref)  : %.3f ms\n",
               time_ms([&]{ spmv(v_cpu, K_ref, u_h); }, reps));

        printf("CPU matrix-free       : %.3f ms\n",
               time_ms([&]{ matvec_cpu_serial(v_cpu, u_h, m); }, reps));

        printf("GPU atomic            : %.3f ms\n",
               time_ms([&]{ matvec_gpu_atomic(d_v, d_u, m.d_coords, m.d_conn, m.nelems); }, reps));

        printf("GPU coloring          : %.3f ms\n",
               time_ms([&]{ matvec_gpu_color(d_v, d_u, m.d_coords, m.d_conn, d_color_groups); }, reps));
    }
    Kokkos::finalize();
    return 0;
}

