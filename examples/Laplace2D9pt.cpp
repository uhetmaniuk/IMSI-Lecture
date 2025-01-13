#include <iostream>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"

using accelerator_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
using host_execution_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;

using real = float;

int main(int argc, char *argv[]) {

    Kokkos::initialize(argc, argv);
    {
        int nx = 32;
        const int ny = nx;
        auto const h = real(1) / real(nx + 1);

        /*
        std::vector<int> mesh, perm;
        std::vector<int> rowind, colptr;
        std::vector<Scalar> nzvals;

        nnodes = nx * ny;

        mesh.resize(nnodes);
        for (i = 0; i < nnodes; i++) {
            mesh[i] = i;
        }

        colptr.resize(nnodes + 1, 0);
        rowind.reserve(9 * nnodes);
        nzvals.reserve(9 * nnodes);

        colptr[0] = 0;
        node = 0;

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {

                if (j > 0) {
                    if (i > 0) {
                        rowind.push_back(mesh(i - 1, j - 1));
                        nzvals.push_back(Scalar(-0.25) / (h * h));
                    }
                    rowind.push_back(mesh(i, j - 1));
                    nzvals.push_back(Scalar(-0.5) / (h * h));
                    if (i + 1 < nx) {
                        rowind.push_back(mesh(i + 1, j - 1));
                        nzvals.push_back(Scalar(-0.25) / (h * h));
                    }
                }

                if (i > 0) {
                    rowind.push_back(mesh(i - 1, j));
                    nzvals.push_back(Scalar(-0.5) / (h * h));
                }
                rowind.push_back(mesh(i, j));
                nzvals.push_back(Scalar(3.0) / (h * h));
                if (i + 1 < nx) {
                    rowind.push_back(mesh(i + 1, j));
                    nzvals.push_back(Scalar(-0.5) / (h * h));
                }

                if (j + 1 < ny) {
                    if (i > 0) {
                        rowind.push_back(mesh(i - 1, j + 1));
                        nzvals.push_back(Scalar(-0.25) / (h * h));
                    }
                    rowind.push_back(mesh(i, j + 1));
                    nzvals.push_back(Scalar(-0.5) / (h * h));
                    if (i + 1 < nx) {
                        rowind.push_back(mesh(i + 1, j + 1));
                        nzvals.push_back(Scalar(-0.25) / (h * h));
                    }
                }

                colptr[node + 1] = int(rowind.size());
                node++;
            }
        }
*/
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
                colIdx(rowPtr(i) + 2) = i + 1;
                values(rowPtr(i) + 2) = -real(1.0 / h);
                rowPtr(i + 1) = rowPtr(i) + 3;
            } else {
                rowPtr(i + 1) = rowPtr(i) + 2;
            }
        }

        Tacho::CrsMatrixBase<real, host_execution_type> h_A;
        h_A.setExternalMatrix(n, n, nnz, rowPtr, colIdx, values);

        Tacho::CrsMatrixBase<real, accelerator_type> A;
        A.createMirror(h_A);

        Tacho::Driver<real, accelerator_type> solver;
        solver.analyze(A.NumRows(), A.RowPtr(), A.Cols());

        // Create numeric tools and levelset tools
        Kokkos::Timer timer;
        solver.initialize();
        double initi_time = timer.seconds();

        //
        // Symbolic structure can be reused
        //
        const int nfacts = 16;
        timer.reset();
        for (int i = 0; i < nfacts; ++i) {
            solver.factorize(A.Values());
        }
        double facto_time = timer.seconds();

        //
        // Tacho::Driver::solve expects a "matrix" for parameters.
        // So even if one linear system with a single RHS needs to be solved,
        // the current interface requires to use a "matrix" with 1 column.
        //
        Kokkos::View<real **, Kokkos::LayoutLeft, typename accelerator_type::execution_space> b("rhs", n, 1);
        Kokkos::View<real **, Kokkos::LayoutLeft, typename accelerator_type::execution_space> x("solution", n, 1);
        Kokkos::View<real **, Kokkos::LayoutLeft, typename accelerator_type::execution_space> wt("workspace", n, 1);

        // Fill the right hand side with random values
        Kokkos::Random_XorShift64_Pool<typename accelerator_type::execution_space> random(13718);
        Kokkos::fill_random(b, random, real(1));

        double solve_time = 0.0;
        const int nsolves = 6;
        for (int i = 0; i < nsolves; ++i) {
            timer.reset();
            solver.solve(x, b, wt);
            solve_time += timer.seconds();
            const double res = solver.computeRelativeResidual(A.Values(), x, b);
            std::cout << "TachoSolver: residual = " << res << "\n";
        }
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << " Initialize Time " << initi_time << std::endl;
        std::cout << " Avg. Facto Time " << facto_time / (double) nfacts << std::endl;
        std::cout << " Avg. Solve Time " << solve_time / (double) nsolves << std::endl;
        std::cout << std::endl;
        solver.release();
    }

    Kokkos::finalize();
    return 0;
}

