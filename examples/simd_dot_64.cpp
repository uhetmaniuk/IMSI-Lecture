// simd_dot_64.cpp
// Kokkos C++17 translation of simd_dot_64.jl
//
// Build example (OpenMP backend):
//   g++ -std=c++17 -O3 -fopenmp -I<kokkos_include> simd_dot_64.cpp \
//       -L<kokkos_lib> -lkokkos -o simd_dot_64
//
// Or with CMake + Kokkos:
//   find_package(Kokkos REQUIRED)
//   target_link_libraries(simd_dot_64 Kokkos::kokkos)

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

// ---------------------------------------------------------------------------
// Naive sequential dot product (mirrors dot_naive in Julia).
// Runs entirely on the host without Kokkos parallelism so the comparison
// remains fair against the @simd-free Julia baseline.
// ---------------------------------------------------------------------------
__attribute__((optimize("O2,no-tree-vectorize")))
double dot_naive(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    const std::size_t n = a.size();

#pragma clang loop vectorize(disable) interleave(disable)
    for (std::size_t i = 0; i < n; ++i) {
        s += a[i] * b[i];
    }
    return s;
}

// ---------------------------------------------------------------------------
// Explicit SIMD dot product using Kokkos::Experimental::native_simd<double>.
//
// This is the closest C++ equivalent of Julia's @simd annotation: instead of
// relying on the compiler to auto-vectorise a scalar loop, we explicitly
// operate on SIMD registers of width simd_type::size().
//
// Algorithm
// ─────────
//   1. Process the array in chunks of simd_type::size() using simd_type loads
//      and fused multiply-add, accumulating into a simd_type register `acc`.
//   2. Reduce `acc` to a scalar with Kokkos::Experimental::reduce().
//   3. Handle the tail (n % simd_width != 0) with a plain scalar loop over
//      the remaining < simd_width elements.
//
// native_simd<double> maps to the widest hardware register available at
// compile time (SSE2: 2 lanes, AVX2: 4 lanes, AVX-512: 8 lanes).
// Compile with -march=native (or the equivalent Kokkos_ARCH_* flag) to get
// the full register width for your CPU.
//
// Note: Kokkos::Experimental::native_simd is a host-only type; this function
// runs on the CPU regardless of the default execution space.
// ---------------------------------------------------------------------------
template<typename Scalar>
Scalar dot_simd(const Scalar* a, const Scalar* b, int n)
{
    using real      = Scalar;
    using simd_type = Kokkos::Experimental::native_simd<real>;
    using tag       = Kokkos::Experimental::element_aligned_tag;

    const int width = static_cast<int>(simd_type::size());
    const int body  = (n / width) * width;  // largest multiple of width ≤ n
    simd_type acc(0.0);
    simd_type va, vb;

    // ── Full chunks: load width doubles at a time, FMA into acc ──────────────
    for (int i = 0; i < body; i += width) {
        va.copy_from(a + i, tag{});
        vb.copy_from(b + i, tag{});
        acc += va * vb;
    }

    // ── Horizontal reductionacross SIMD lanes ────────────────────────────────
    real result = 0.0;
    for (int k = 0; k < width; ++k) {
        result += acc[k];
    }

    // ── Tail: plain scalar loop for the remaining < width elements ────────────
    for (int i = body; i < n; ++i)
        result += a[i] * b[i];

    return result;
}

// Convenience overload for std::vector
double dot_simd(const std::vector<double>& a, const std::vector<double>& b)
{
    return dot_simd(a.data(), b.data(), static_cast<int>(a.size()));
}

// ---------------------------------------------------------------------------
// Kokkos parallel_reduce dot product (mirrors dot_simd / @simd in Julia).
//
// On a CPU backend (Serial, OpenMP, Threads) Kokkos will emit SIMD-friendly
// loops.  On a GPU backend (CUDA, HIP) it maps to a standard reduction
// kernel.  The ExecutionSpace template parameter lets you switch backends at
// compile time without touching the algorithm.
// ---------------------------------------------------------------------------
template <class ExecutionSpace>
double dot_kokkos(const Kokkos::View<const double*, ExecutionSpace>& a,
                  const Kokkos::View<const double*, ExecutionSpace>& b)
{
    const int n = static_cast<int>(a.extent(0));
    double result = 0.0;

    Kokkos::parallel_reduce(
        "dot_kokkos",
        Kokkos::RangePolicy<ExecutionSpace>(0, n),
        KOKKOS_LAMBDA(const int i, double& sum) {
            sum += a(i) * b(i);
        },
        result
    );

    return result;
}

// ---------------------------------------------------------------------------
// Simple timer wrapper around Kokkos::Timer.
// Returns elapsed seconds for a single functor call.
// ---------------------------------------------------------------------------
template <class F>
double time_seconds(F&& fn, int repeats = 5)
{
    // warm-up
    fn();
    Kokkos::fence();

    Kokkos::Timer timer;
    for (int r = 0; r < repeats; ++r) {
        fn();
    }
    Kokkos::fence();
    return timer.seconds() / repeats;
}

// ---------------------------------------------------------------------------
// Helpers to pretty-print timings (mimics BenchmarkTools @btime output)
// ---------------------------------------------------------------------------
void print_time(const char* label, double seconds)
{
    if (seconds < 1e-6)
        std::printf("  %-26s  %.3f ns\n", label, seconds * 1e9);
    else if (seconds < 1e-3)
        std::printf("  %-26s  %.3f μs\n", label, seconds * 1e6);
    else if (seconds < 1.0)
        std::printf("  %-26s  %.3f ms\n", label, seconds * 1e3);
    else
        std::printf("  %-26s  %.3f s\n",  label, seconds);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
        // Generate the same log-spaced size sequence as Julia:
        //   p = 10 .^ (2:0.1:6)  →  n in [100, 1_000_000]
        std::vector<int> sizes;
        for (double exp = 2.0; exp <= 6.0 + 1e-9; exp += 0.1) {
            int n = static_cast<int>(std::ceil(std::pow(10.0, exp)));
            sizes.push_back(n);
        }

        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        std::printf("%-6s  %-12s  %-14s  %-14s\n",
                    "n", "rel_err_naive", "rel_err_kokkos", "");
        std::printf("%s\n", std::string(60, '-').c_str());

        for (int n : sizes) {
            // -----------------------------------------------------------------
            // Fill host vectors with random data (mirrors rand(n) in Julia)
            // -----------------------------------------------------------------
            std::vector<double> h_a(n), h_b(n);
            for (int i = 0; i < n; ++i) {
                h_a[i] = dist(rng);
                h_b[i] = dist(rng);
            }

            // -----------------------------------------------------------------
            // Reference: std::inner_product (high-accuracy, mirrors BLAS dot)
            // -----------------------------------------------------------------
            double s_ref = 0.0;
            for (int i = 0; i < n; ++i) s_ref += h_a[i] * h_b[i];

            // -----------------------------------------------------------------
            // Naive dot (host, no SIMD)
            // -----------------------------------------------------------------
            double s_naive = dot_naive(h_a, h_b);

            // -----------------------------------------------------------------
            // Explicit SIMD dot (host, native_simd<double>)
            // -----------------------------------------------------------------
            double s_simd = dot_simd(h_a, h_b);

            // -----------------------------------------------------------------
            // Copy data to Kokkos Views (device memory when using GPU backend)
            // -----------------------------------------------------------------
            using DevView = Kokkos::View<double*, Kokkos::DefaultExecutionSpace>;
            using ConstDevView = Kokkos::View<const double*, Kokkos::DefaultExecutionSpace>;

            DevView d_a("a", n), d_b("b", n);
            auto h_a_view = Kokkos::create_mirror_view(d_a);
            auto h_b_view = Kokkos::create_mirror_view(d_b);
            for (int i = 0; i < n; ++i) {
                h_a_view(i) = h_a[i];
                h_b_view(i) = h_b[i];
            }
            Kokkos::deep_copy(d_a, h_a_view);
            Kokkos::deep_copy(d_b, h_b_view);
            ConstDevView ca = d_a, cb = d_b;

            double s_kokkos = dot_kokkos<Kokkos::Serial>(ca, cb);

            // -----------------------------------------------------------------
            // Relative errors (mirrors the @sprintf line in Julia)
            // -----------------------------------------------------------------
            double rel_naive  = (s_ref != 0.0) ? (s_naive  - s_ref) / s_ref : 0.0;
            double rel_simd   = (s_ref != 0.0) ? (s_simd   - s_ref) / s_ref : 0.0;
            double rel_kokkos = (s_ref != 0.0) ? (s_kokkos - s_ref) / s_ref : 0.0;
            std::printf("n=%-10d  naive: %+.6e  simd: %+.6e  kokkos: %+.6e\n",
                        n, rel_naive, rel_simd, rel_kokkos);

            // -----------------------------------------------------------------
            // Benchmarks (mirrors @btime in Julia)
            // -----------------------------------------------------------------
            const int repeats = 1024;

            double t_naive = time_seconds([&]{ volatile double v = dot_naive(h_a, h_b); (void)v; }, repeats);
            print_time("Naive (host scalar):", t_naive);

            double t_simd = time_seconds([&]{ volatile double v = dot_simd(h_a, h_b); (void)v; }, repeats);
            print_time("SIMD (native_simd):", t_simd);

            double t_kokkos = time_seconds([&]{ volatile double v = dot_kokkos<Kokkos::Serial>(ca, cb); (void)v; }, repeats);
            print_time("Kokkos parallel_reduce:", t_kokkos);

            std::printf("\n");
        }
    }
    Kokkos::finalize();
    return 0;
}

