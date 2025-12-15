#include <immintrin.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_SIMD.hpp>
#include <cmath>
#include <iostream>
#include <vector>

int
main(int argc, char** argv)
{
  const int rep = 8192;
  double    sum = 0.;
  using real    = float;

  Kokkos::initialize(argc, argv);
  {
    using simd_type = Kokkos::Experimental::native_simd<real>;
    std::cout << "double Kokkos " << Kokkos::Experimental::native_simd<double>::size() << "\n";
    std::cout << " float Kokkos " << Kokkos::Experimental::native_simd<float>::size() << "\n";
    std::cout << " int32 Kokkos " << Kokkos::Experimental::native_simd<int32_t>::size() << "\n";
    std::cout << " int64 Kokkos " << Kokkos::Experimental::native_simd<int64_t>::size() << "\n";

    std::cout << "\n";
    std::cout << " p   |   for Loop   for Comp.      Kokkos";
    std::cout << "\n";

    for (int pk = 32; pk < 80000; pk *= 2) {
      for (int p = pk - 2; p <= pk + 2; p += 1) {
        Kokkos::Random_XorShift64_Pool<> random_pool(/*seed=*/12345);
        auto                             generator = random_pool.get_state();

        std::vector<real> x(p), y(p), z(p);
        for (int i = 0; i < p; ++i) {
          x[i] = generator.drand(0., 1.);
          y[i] = generator.drand(0., 1.);
          z[i] = generator.drand(0., 1.);
        }

        std::cout << p << " | ";

        std::vector<real> r(p);
        sum = 0.0;
        std::chrono::duration<double> dt_ref(0);
        for (int ir = 0; ir < rep; ++ir) {
          auto start = std::chrono::high_resolution_clock::now();
#pragma clang loop vectorize(disable) interleave(disable)
          for (int jj = 0; jj < p; ++jj) {
            r[jj] = x[jj] * x[jj] + y[jj] * y[jj] + z[jj] * z[jj];
          }
          auto end = std::chrono::high_resolution_clock::now();
          dt_ref += end - start;
          sum += r[0] + r[p - 1];
        }
        std::cout << "      1      ";
        //--- Use for debugging
        // std::cout << " sum " << sum << "\n";

        //
        // LLVM vectorization https://llvm.org/docs/Vectorizers.html
        //
        sum = 0.0;
        r.assign(p, 0.0);
        std::chrono::duration<double> dt88(0);
        for (int ir = 0; ir < rep; ++ir) {
          auto start = std::chrono::high_resolution_clock::now();
#pragma clang loop vectorize(enable) interleave(enable)
          for (int jj = 0; jj < p; ++jj) {
            r[jj] = x[jj] * x[jj] + y[jj] * y[jj] + z[jj] * z[jj];
          }
          auto end = std::chrono::high_resolution_clock::now();
          dt88 += end - start;
          sum += r[0] + r[p - 1];
        }
        std::cout << dt_ref.count() / dt88.count() << "   ";
        //--- Use for debugging
        // std::cout << " sum " << sum << "\n";

        /*
                        sum = 0.0;
                        r.assign(r.size(), 0.0);
                        std::chrono::duration<double> dt44(0);
                        for (int ir = 0; ir < rep; ++ir) {
                            auto start = std::chrono::high_resolution_clock::now();
                            int i = 0;
                            if constexpr (std::is_same_v<real, double>) {
                                for (; i + Vc::double_v::size() < p; i += Vc::double_v::size()) {
                                    Vc::double_v v_x(&x[i]), v_y(&y[i]), v_z(&z[i]);
                                    auto v_r = v_x * v_x + v_y * v_y + v_z * v_z;
                                    v_r.store(&r[i]);
                                }
                            } else if constexpr (std::is_same_v<real, float>) {
                                for (; i + Vc::float_v::size() < p; i += Vc::float_v::size()) {
                                    Vc::float_v v_x(&x[i]), v_y(&y[i]), v_z(&z[i]);
                                    auto v_r = v_x * v_x + v_y * v_y + v_z * v_z;
                                    v_r.store(&r[i]);
                                }
                            }
                            for (; i < p; ++i) {
                                r[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
                            }
                            auto end = std::chrono::high_resolution_clock::now();
                            dt44 += end - start;
                            sum += r[0] + r[p - 1];
                        }
                        std::cout << "  " << dt_ref.count() / dt44.count() << "  ";
                        //--- Use for debugging
                        // std::cout << " sum " << sum << "\n";
        */

        using tag_type      = Kokkos::Experimental::element_aligned_tag;
        constexpr int width = int(simd_type::size());
        sum                 = 0.0;
        r.assign(p, 0.0);
        std::chrono::duration<double> dt33(0);
        for (int ir = 0; ir < rep; ++ir) {
          auto      start = std::chrono::high_resolution_clock::now();
          simd_type sx, sy, sz;
          int       i = 0;
          for (; i + width <= p; i += width) {
            sx.copy_from(&x[i], Kokkos::Experimental::element_aligned_tag());
            sy.copy_from(&y[i], Kokkos::Experimental::element_aligned_tag());
            sz.copy_from(&z[i], Kokkos::Experimental::element_aligned_tag());
            auto sr = sx * sx + sy * sy + sz * sz;
            sr.copy_to(&r[i], Kokkos::Experimental::element_aligned_tag());
          }
          for (; i < p; ++i) {
            r[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
          }
          auto end = std::chrono::high_resolution_clock::now();
          dt33 += end - start;
          sum += r[0] + r[p - 1];
        }
        std::cout << dt_ref.count() / dt33.count() << "  ";
        //--- Use for debugging
        // std::cout << " sum " << sum << "\n";

        sum = 0.0;
        r.assign(p, 0.0);
        auto*                         rPtr = r.data();
        auto*                         xPtr = x.data();
        auto*                         yPtr = y.data();
        auto*                         zPtr = z.data();
        std::chrono::duration<double> dtSerial(0);
        for (int ir = 0; ir < rep; ++ir) {
          auto start = std::chrono::high_resolution_clock::now();
          Kokkos::parallel_for(
              "Serial", Kokkos::RangePolicy<Kokkos::Serial>(0, p), [=](const int k) {
                rPtr[k] = xPtr[k] * xPtr[k] + yPtr[k] * yPtr[k] + zPtr[k] * zPtr[k];
              });
          auto end = std::chrono::high_resolution_clock::now();
          dtSerial += end - start;
          sum += r[0] + r[p - 1];
        }
        std::cout << dt_ref.count() / dtSerial.count() << "  ";
        //--- Use for debugging
        // std::cout << " sum " << sum << "\n";

        sum = 0.0;
        r.assign(p, 0.0);
        std::chrono::duration<double> dtOpenMP(0);
        for (int ir = 0; ir < rep; ++ir) {
          auto start = std::chrono::high_resolution_clock::now();
          Kokkos::parallel_for("OMP", Kokkos::RangePolicy<Kokkos::OpenMP>(0, p), [=](const int k) {
            rPtr[k] = xPtr[k] * xPtr[k] + yPtr[k] * yPtr[k] + zPtr[k] * zPtr[k];
          });
          auto end = std::chrono::high_resolution_clock::now();
          dtOpenMP += end - start;
          sum += r[0] + r[p - 1];
        }
        std::cout << dt_ref.count() / dtOpenMP.count() << "  ";
        //--- Use for debugging
        // std::cout << " sum " << sum << "\n";

        if (p % 4 == 0) {
          sum = 0.0;
          r.assign(p, 0.0);
          std::chrono::duration<double> dt44(0);
          // TeamPolicy(League Size, Team Size, Vector Length)
          Kokkos::TeamPolicy<Kokkos::Serial> policy(p / 4, 1, 4);
          for (int ir = 0; ir < rep; ++ir) {
            auto start = std::chrono::high_resolution_clock::now();
            Kokkos::parallel_for(
                policy, [&](const Kokkos::TeamPolicy<Kokkos::Serial>::member_type& team_member) {
                  int k = team_member.league_rank();
                  Kokkos::parallel_for(
                      Kokkos::ThreadVectorRange(team_member, 4 * k, 4 * k + 4), [&](const int j) {
                        r[j] = x[j] * x[j] + y[j] * y[j] + z[j] * z[j];
                      });
                });
            auto end = std::chrono::high_resolution_clock::now();
            dt44 += end - start;
            sum += r[0] + r[p - 1];
          }
          std::cout << dt_ref.count() / dt44.count() << "  ";
        }

        if (p % 4 == 0) {
          sum = 0.0;
          r.assign(p, 0.0);
          std::chrono::duration<double> dt44(0);
          // TeamPolicy(League Size, Team Size, Vector Length)
          Kokkos::TeamPolicy<Kokkos::OpenMP> policy(p / 4, 1, 4);
          for (int ir = 0; ir < rep; ++ir) {
            auto start = std::chrono::high_resolution_clock::now();
            Kokkos::parallel_for(
                policy, [&](const Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type& team_member) {
                  int k = team_member.league_rank();
                  Kokkos::parallel_for(
                      Kokkos::ThreadVectorRange(team_member, 4 * k, 4 * k + 4), [&](const int j) {
                        r[j] = x[j] * x[j] + y[j] * y[j] + z[j] * z[j];
                      });
                });
            auto end = std::chrono::high_resolution_clock::now();
            dt44 += end - start;
            sum += r[0] + r[p - 1];
          }
          std::cout << dt_ref.count() / dt44.count() << "  ";
        }

        std::cout << std::endl;
      }
    }
  }
  Kokkos::finalize();
}
