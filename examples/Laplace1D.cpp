#include <iostream>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"

int main(int argc, char *argv[]) {

    Kokkos::initialize(argc, argv);
    {
        using execution_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
        using host_execution_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;

        using real = double; // float;

        int n = 64;
        int nnz = 3 * n - 2;

        real h = real(1) / n;

        Kokkos::View<real *> values("values", nnz);
        Kokkos::View<size_t *> rowPtr("row pointer", n + 1);
        Kokkos::View<int *> colIdx("column indices", nnz);

        rowPtr(0) = 0;
        rowPtr(1) = 2;
        colIdx(0) = 0;
        colIdx(1) = 1;
        values(0) = real(2.0 / h);
        values(1) = -real(1.0 / h);
        for (int i = 1; i < n; ++i) {
            colIdx(rowPtr(i)) = i - 1;
            colIdx(rowPtr(i) + 1) = i;
            values(rowPtr(i)) = -real(1.0 / h);
            values(rowPtr(i) + 1) = real(2.0 / h);
            if (i + 1 < n) {
                rowPtr(i + 1) = rowPtr(i) + 3;
                colIdx(rowPtr(i) + 2) = i + 1;
                values(rowPtr(i) + 2) = -real(1.0 / h);
            } else {
                rowPtr(i + 1) = rowPtr(i) + 2;
            }
        }

        Tacho::CrsMatrixBase<real, execution_type> A;
        A.setExternalMatrix(n, n, nnz, rowPtr, colIdx, values);;

        Tacho::Driver<real, execution_type> solver;

        auto values_on_executor = Kokkos::create_mirror_view(typename execution_type::memory_space(), A.Values());
        Kokkos::deep_copy(values_on_executor, A.Values());

        solver.analyze(A.NumRows(), A.RowPtr(), A.Cols());

        /// create numeric tools and levelset tools
        Kokkos::Timer timer;
        solver.initialize();
        double initi_time = timer.seconds();

        /// symbolic structure can be reused
        const int nfacts = 32;
        timer.reset();
        for (int i = 0; i < nfacts; ++i) {
            solver.factorize(values_on_executor);
        }
        double facto_time = timer.seconds();

        Kokkos::View<real **, Kokkos::LayoutLeft, typename execution_type::execution_space> b("rhs", n, 1);
        Kokkos::View<real **, Kokkos::LayoutLeft, typename execution_type::execution_space> x("solution", n, 1);
        Kokkos::View<real **, Kokkos::LayoutLeft, typename execution_type::execution_space> wt("workspace", n, 1);

        Kokkos::Random_XorShift64_Pool<typename execution_type::execution_space> random(13718);
        Kokkos::fill_random(b, random, real(1));

        double solve_time = 0.0;
        const int nsolves = 6;
        for (int i = 0; i < nsolves; ++i) {
            timer.reset();
            solver.solve(x, b, wt);
            solve_time += timer.seconds();
            const double res = solver.computeRelativeResidual(values_on_executor, x, b);
            std::cout << "TachoSolver: residual = " << res << "\n";
        }
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << " Initi Time " << initi_time << std::endl;
        std::cout << " Avg. Facto Time " << facto_time / (double) nfacts << std::endl;
        std::cout << " Solve Time " << solve_time / (double) nsolves << std::endl;
        std::cout << std::endl;
        solver.release();
    }

    Kokkos::finalize();
    return 0;
}

