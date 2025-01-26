#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_SIMD.hpp>

#include <Vc/Vc>

#include "../TPL/eigen/Eigen/Dense"

#include <cmath>
#include <immintrin.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {

    const int nNode = 16;
    const int dim = 2;
    const int rep = 4096;
    double sum = 0.;

    Kokkos::initialize(argc,argv);
    {
        using simd_type = Kokkos::Experimental::native_simd<double>;
        std::cout << "double " << Kokkos::Experimental::native_simd<double>::size() << "\n";
        std::cout << "float " << Kokkos::Experimental::native_simd<float>::size() << "\n";
        std::cout << "int32 " << Kokkos::Experimental::native_simd<int32_t>::size() << "\n";
        std::cout << "int64 " << Kokkos::Experimental::native_simd<int64_t>::size() << "\n";
        //
        std::vector<double> g(2 * nNode);
        Eigen::Matrix<double, nNode, nNode, Eigen::ColMajor> kele;
        auto start = std::chrono::high_resolution_clock::now();
        for (int ir = 0; ir < rep; ++ir) {
          for (int in = 0; in < 2; ++in) {
            for (int jn = 0; jn < nNode; ++jn) {
              g[jn + in * nNode] = in + jn + ir;
            }
          }
          double ax = 4.0, ay = 5.0;
          for (int j = 0; j < nNode; j += 1) {
            simd_type a(ax);
            simd_type gj(g[j]);
            for (int i = 0; i < nNode; i += Kokkos::Experimental::native_simd<double>::size()) {
              simd_type col;
              col.copy_from(&g[i], Kokkos::Experimental::simd_flag_default);
              col *= a*gj;
              col.copy_to(&kele(i, j), Kokkos::Experimental::simd_flag_default);
            }
          }
          //
          for (int j = 0; j < nNode; ++j) {
            simd_type a(ay);
            simd_type hj(g[j + nNode]);
            for (int i = 0; i < nNode; i += Kokkos::Experimental::native_simd<double>::size()) {
              simd_type dy;
              dy.copy_from(&g[i + nNode], Kokkos::Experimental::simd_flag_default);
              simd_type col;
              col.copy_from(&kele(i, j), Kokkos::Experimental::simd_flag_default);
              col += dy * a * hj;
              col.copy_to(&kele(i, j), Kokkos::Experimental::simd_flag_default);
            }
          }
          sum += kele(0, nNode - 2) + kele(nNode - 1, nNode - 1);
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = end - start;
        std::cout << " Kokkos loop " << dt.count() << "\n";
        std::cout << " sum " << sum << "\n";

std::cout << " Vc float " << Vc::float_v::size()
<< " double " << Vc::double_v::size()
<< "\n";
using real = double;

for (int p = 512; p < 3000; p *= 2) {
    Kokkos::Random_XorShift64_Pool<> random_pool(/*seed=*/12345);
    auto generator = random_pool.get_state();

    std::vector<real> x(p), y(p), z(p);
    for (int i = 0; i < p; ++i) {
      x[i] = generator.drand(0., 1.);
      y[i] = generator.drand(0., 1.);
      z[i] = generator.drand(0., 1.);
    }

    std::vector<real> r(p);
    sum = 0.0;
    std::chrono::duration<double> dt22(0);
    for (int ir = 0; ir < rep; ++ir) {
    start = std::chrono::high_resolution_clock::now();
#pragma clang loop vectorize(disable) interleave(disable)
      for (int jj = 0; jj < p; ++jj) {
         r[jj] = sqrt(x[jj] * x[jj] + sin(y[jj]) * sin(y[jj]) + exp(z[jj]) * z[jj]);
         //r[jj] = x[jj] * x[jj] + y[jj] * y[jj] + z[jj] * z[jj];
      }
    end = std::chrono::high_resolution_clock::now();
    dt22 += end - start;
      sum += r[0] + r[p-1];
    }
    std::cout << p << " Manual loop " << dt22.count() / (p * rep) << "\n";
    std::cout << " sum " << sum << "\n";

    sum = 0.0;
    std::chrono::duration<double> dt88(0);
    for (int ir = 0; ir < rep; ++ir) {
        start = std::chrono::high_resolution_clock::now();
#pragma clang loop vectorize(enable) interleave(enable)
        for (int jj = 0; jj < p; ++jj) {
            r[jj] = sqrt(x[jj] * x[jj] + sin(y[jj]) * sin(y[jj]) + exp(z[jj]) * z[jj]);
            //r[jj] = x[jj] * x[jj] + y[jj] * y[jj] + z[jj] * z[jj];
        }
        end = std::chrono::high_resolution_clock::now();
        dt88 += end - start;
        sum += r[0] + r[p-1];
    }
    std::cout << p << " Manual loop " << dt88.count() / (p * rep) << "\n";
    std::cout << " sum " << sum << "\n";

    sum = 0.0;
    r.assign(r.size(), 0.0);
    std::chrono::duration<double> dt44(0);
    for (int ir = 0; ir < rep; ++ir) {
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < p; i += Vc::double_v::size()) {
      Vc::double_v v_x(&x[i]), v_y(&y[i]), v_z(&z[i]);
      auto v_r = sqrt(v_x * v_x + Vc::sin(v_y) * Vc::sin(v_y) + Vc::exp(v_z) * v_z);
      //auto v_r = v_x * v_x + v_y * v_y + v_z * v_z;
      v_r.store(&r[i]);
    }
    end = std::chrono::high_resolution_clock::now();
    dt44 += end - start;
      sum += r[0] + r[p-1];
    }
    std::cout << p << " Vc loop " << dt44.count() / (p * rep) << "\n";
    std::cout << " sum " << sum << "\n";

using tag_type = Kokkos::Experimental::element_aligned_tag;
constexpr int width = int(simd_type::size());
    sum = 0.0;
    std::chrono::duration<double> dt33(0);
    for (int ir = 0; ir < rep; ++ir) {
    start = std::chrono::high_resolution_clock::now();
  simd_type sx, sy, sz;
for (int i = 0; i < p; i += simd_type::size()) {
  sx.copy_from(&x[i], Kokkos::Experimental::element_aligned_tag());
  sy.copy_from(&y[i], Kokkos::Experimental::element_aligned_tag());
  sz.copy_from(&z[i], Kokkos::Experimental::element_aligned_tag());
  auto sr = Kokkos::sqrt(sx * sx + Kokkos::sin(sy) * Kokkos::sin(sy) + Kokkos::exp(sz) * sz);
  //auto sr = sx * sx + sy * sy + sz * sz;
  sr.copy_to(&r[i], Kokkos::Experimental::element_aligned_tag());
}
    end = std::chrono::high_resolution_clock::now();
    dt33 += end - start;
sum += r[0] + r[p-1];
}
    std::cout << " Kokkos SIMD " << dt33.count() / (p * rep) << "\n";
    std::cout << " sum " << sum << "\n";

}


    }
    Kokkos::finalize();

}