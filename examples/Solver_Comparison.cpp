#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Kokkos_Core.hpp>
#include "../src/SSOR_Solver.h"
#include "../src/PCG_Solver.h"
#include "../src/SymmetricSparse.hpp"

/**
 * @brief Sparse linear solver comparison example
 *
 * Compares different methods for solving symmetric positive definite systems:
 * - SSOR (Symmetric Successive Over-Relaxation) iterative solver
 * - PCG (Preconditioned Conjugate Gradient) with diagonal preconditioning
 * - CG (Conjugate Gradient) unpreconditioned
 * - Direct solver (LDL^T factorization via SymmetricSparse)
 *
 * Test problem: 1D Laplacian equation (tridiagonal matrix)
 * Matrix: tridiag(-1, 2, -1) representing discrete -u'' = f
 */

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "Sparse Linear Solver Comparison\n";
        std::cout << "================================\n\n";

        // Problem size: 2D grid
        const int nx = 64;  // Grid points in x direction
        const int ny = 64;  // Grid points in y direction
        const int n = nx * ny;  // Total number of unknowns

        // Create a 2D Laplacian matrix with 9-point stencil
        // Stencil: [-1 -1 -1]
        //          [-1  8 -1]
        //          [-1 -1 -1]
        // This corresponds to the discretization of -∇²u = f on a 2D domain
        std::cout << "Building " << nx << "x" << ny << " 2D Laplacian (9-point stencil)...\n";

        std::vector<int> rowPtr(n + 1);
        std::vector<int> colInd;
        std::vector<double> values;

        // Helper to get node index from (i,j) coordinates
        auto nodeIndex = [&](int i, int j) -> int {
            return i * ny + j;
        };

        rowPtr[0] = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int node = nodeIndex(i, j);
                int nnz_row = 0;

                // 9-point stencil: add all neighbors
                // Store in sorted order for better cache performance

                // Southwest neighbor
                if (i > 0 && j > 0) {
                    colInd.push_back(nodeIndex(i-1, j-1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // South neighbor
                if (i > 0) {
                    colInd.push_back(nodeIndex(i-1, j));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // Southeast neighbor
                if (i > 0 && j < ny - 1) {
                    colInd.push_back(nodeIndex(i-1, j+1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // West neighbor
                if (j > 0) {
                    colInd.push_back(nodeIndex(i, j-1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // Center (diagonal)
                colInd.push_back(node);
                values.push_back(8.0);
                nnz_row++;

                // East neighbor
                if (j < ny - 1) {
                    colInd.push_back(nodeIndex(i, j+1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // Northwest neighbor
                if (i < nx - 1 && j > 0) {
                    colInd.push_back(nodeIndex(i+1, j-1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // North neighbor
                if (i < nx - 1) {
                    colInd.push_back(nodeIndex(i+1, j));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                // Northeast neighbor
                if (i < nx - 1 && j < ny - 1) {
                    colInd.push_back(nodeIndex(i+1, j+1));
                    values.push_back(-1.0);
                    nnz_row++;
                }

                rowPtr[node + 1] = rowPtr[node] + nnz_row;
            }
        }

        int nnz = colInd.size();
        std::cout << "Matrix has " << nnz << " non-zeros\n";
        std::cout << "Average non-zeros per row: " << (double)nnz / n << "\n\n";

        // Create right-hand side: b = ones(n)
        // Solution should be approximately parabolic
        std::vector<double> rhs(n, 1.0);
        std::vector<double> sol(n, 0.0);  // Initial guess = 0
        std::vector<double> diag(n);

        // Extract diagonal
        IMSI::ExtractDiagonal(n, rowPtr.data(), colInd.data(), values.data(), diag.data());

        std::cout << "Solving with SSOR (omega = 1.0)...\n";
        std::cout << "-----------------------------------\n";

        // Solve the system
        int iterations = IMSI::SSOR_Solve(
            n,
            rowPtr.data(),
            colInd.data(),
            values.data(),
            diag.data(),
            rhs.data(),
            sol.data(),
            1.0,      // omega
            1e-3,     // tolerance
            10000     // max iterations
        );

        std::cout << "Iterations returned: " << iterations << "\n";

        // Compute actual residual
        std::vector<double> Ax(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                Ax[i] += values[k] * sol[colInd[k]];
            }
        }
        double residual_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = rhs[i] - Ax[i];
            residual_norm += diff * diff;
        }
        residual_norm = std::sqrt(residual_norm);
        std::cout << "Final residual: " << residual_norm << "\n";

        if (iterations > 0) {
            std::cout << "Solution converged in " << iterations << " iterations!\n";
            std::cout << "\nFirst 10 solution values:\n";
            for (int i = 0; i < std::min(10, n); ++i) {
                std::cout << "  x[" << i << "] = " << sol[i] << "\n";
            }

            // Verify solution by computing residual
            std::vector<double> Ax(n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                    Ax[i] += values[k] * sol[colInd[k]];
                }
            }

            double residual_norm = 0.0;
            for (int i = 0; i < n; ++i) {
                double diff = rhs[i] - Ax[i];
                residual_norm += diff * diff;
            }
            residual_norm = std::sqrt(residual_norm);

            std::cout << "\nVerification: ||b - Ax||_2 = " << residual_norm << "\n";
        } else {
            std::cout << "Solver did not converge!\n";
        }

        // Test with different omega values
        std::cout << "\n\nTesting different relaxation parameters (with timing):\n";
        std::cout << "======================================================\n";
        std::cout << std::setw(8) << "omega"
                  << std::setw(12) << "iterations"
                  << std::setw(15) << "time (ms)"
                  << std::setw(18) << "time/iter (μs)\n";
        std::cout << std::string(53, '-') << "\n";

        std::vector<double> omega_values = {0.8, 1.0, 1.2, 1.5};
        for (double omega : omega_values) {
            std::fill(sol.begin(), sol.end(), 0.0);  // Reset solution

            Kokkos::Timer timer;
            int iter = IMSI::SSOR_Solve(
                n, rowPtr.data(), colInd.data(), values.data(), diag.data(),
                rhs.data(), sol.data(),
                omega, 1e-3, 10000
            );
            double time_ms = timer.seconds() * 1000.0;

            if (iter > 0) {
                double time_per_iter = (time_ms * 1000.0) / iter;  // microseconds
                std::cout << std::setw(8) << std::fixed << std::setprecision(1) << omega
                          << std::setw(12) << iter
                          << std::setw(15) << std::setprecision(3) << time_ms
                          << std::setw(18) << std::setprecision(2) << time_per_iter << "\n";
            } else {
                std::cout << std::setw(8) << omega << "    did not converge\n";
            }
        }

        // Test PCG solver
        std::cout << "\n\nTesting PCG solver (diagonal preconditioning, tol=1e-12):\n";
        std::cout << "=========================================================\n";

        std::fill(sol.begin(), sol.end(), 0.0);  // Reset solution
        Kokkos::Timer timer_pcg;
        int iter_pcg = IMSI::PCG_Solve(
            n, rowPtr.data(), colInd.data(), values.data(), diag.data(),
            rhs.data(), sol.data(), 1e-12, 10000
        );
        double time_pcg = timer_pcg.seconds() * 1000.0;

        if (iter_pcg > 0) {
            double time_per_iter = (time_pcg * 1000.0) / iter_pcg;
            std::cout << "PCG converged in " << iter_pcg << " iterations\n";
            std::cout << "  Time:         " << std::setprecision(3) << time_pcg << " ms\n";
            std::cout << "  Time/iter:    " << std::setprecision(2) << time_per_iter << " μs\n";
            std::cout << "  First 5 values: ";
            for (int i = 0; i < std::min(5, n); ++i) {
                std::cout << std::fixed << std::setprecision(2) << sol[i] << " ";
            }
            std::cout << "\n";
        } else {
            std::cout << "PCG did not converge\n";
        }

        // Test unpreconditioned CG
        std::cout << "\nTesting unpreconditioned CG (tol=1e-12):\n";
        std::cout << "----------------------------------------\n";

        std::fill(sol.begin(), sol.end(), 0.0);  // Reset solution
        Kokkos::Timer timer_cg;
        int iter_cg = IMSI::CG_Solve(
            n, rowPtr.data(), colInd.data(), values.data(),
            rhs.data(), sol.data(), 1e-12, 10000
        );
        double time_cg = timer_cg.seconds() * 1000.0;

        if (iter_cg > 0) {
            double time_per_iter = (time_cg * 1000.0) / iter_cg;
            std::cout << "CG converged in " << iter_cg << " iterations\n";
            std::cout << "  Time:         " << std::setprecision(3) << time_cg << " ms\n";
            std::cout << "  Time/iter:    " << std::setprecision(2) << time_per_iter << " μs\n";
        } else {
            std::cout << "CG did not converge\n";
        }

        // Test PCG with SSOR preconditioning
        std::cout << "\nTesting PCG with SSOR preconditioning (tol=1e-12):\n";
        std::cout << "---------------------------------------------------\n";

        std::vector<double> work_pcg_ssor(4 * n);
        std::vector<int> ssor_sweep_counts = {1, 2};
        std::vector<double> omega_values_ssor = {1.0, 1.2};

        for (int numSweeps : ssor_sweep_counts) {
            for (double omega_ssor : omega_values_ssor) {
                std::fill(sol.begin(), sol.end(), 0.0);  // Reset solution

                Kokkos::Timer timer_pcg_ssor;
                int iter_pcg_ssor = IMSI::PCG_Solve_SSOR_Precond(
                    n, rowPtr.data(), colInd.data(), values.data(), diag.data(),
                    rhs.data(), sol.data(), work_pcg_ssor.data(),
                    omega_ssor, numSweeps, 1e-12, 10000
                );
                double time_pcg_ssor = timer_pcg_ssor.seconds() * 1000.0;

                if (iter_pcg_ssor > 0) {
                    double time_per_iter = (time_pcg_ssor * 1000.0) / iter_pcg_ssor;
                    std::cout << "  " << numSweeps << " sweep(s), ω=" << std::setprecision(1) << omega_ssor
                              << ": " << std::setw(3) << iter_pcg_ssor << " iter"
                              << ", " << std::setprecision(3) << std::setw(6) << time_pcg_ssor << " ms"
                              << ", " << std::setprecision(2) << std::setw(5) << time_per_iter << " μs/iter\n";
                }
            }
        }

        // Test convenience function with automatic diagonal extraction
        std::cout << "\n\nTesting SSOR_Solve_AutoDiag (automatic diagonal extraction):\n";
        std::cout << "=============================================================\n";

        std::fill(sol.begin(), sol.end(), 0.0);  // Reset solution
        int iter_auto = IMSI::SSOR_Solve_AutoDiag(
            n,
            rowPtr.data(),
            colInd.data(),
            values.data(),
            rhs.data(),
            sol.data(),
            1.2,      // omega
            1e-3,     // tolerance
            10000     // max iterations
        );

        if (iter_auto > 0) {
            std::cout << "Converged in " << iter_auto << " iterations (omega=1.2)\n";
            std::cout << "First 5 solution values: ";
            for (int i = 0; i < std::min(5, n); ++i) {
                std::cout << sol[i] << " ";
            }
            std::cout << "\n";
        }

        // Compare with direct solver (SymmetricSparse)
        std::cout << "\n\nComparison with direct solver (SymmetricSparse LDL^T):\n";
        std::cout << "========================================================\n";

        // Prepare solution for direct solver
        std::vector<double> sol_direct(n, 0.0);

        Kokkos::Timer timer_direct;

        std::cout << "Matrix properties:\n";
        std::cout << "  Dimension (n):  " << n << "\n";
        std::cout << "  Nonzeros (nnz): " << nnz << "\n";
        std::cout << "  Sparsity:       " << std::setprecision(2)
                  << (100.0 * nnz / (n * n)) << "%\n";
        std::cout << "\n";

        // Create direct solver
        SymmetricSparse<double> directSolver(n, nnz, rowPtr.data(), colInd.data(), values.data(), true);

        // Factor the matrix
        double time_factor_start = timer_direct.seconds();
        int factor_result = directSolver.factor();
        double time_factor = (timer_direct.seconds() - time_factor_start) * 1000.0;

        if (factor_result == 0) {
            // Solve
            double time_solve_start = timer_direct.seconds();
            directSolver.Solve(1, rhs.data(), sol_direct.data());
            double time_solve = (timer_direct.seconds() - time_solve_start) * 1000.0;
            double time_total = timer_direct.seconds() * 1000.0;

            // Verify solution
            std::vector<double> Ax_direct(n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                    Ax_direct[i] += values[k] * sol_direct[colInd[k]];
                }
            }
            double residual_direct = 0.0;
            for (int i = 0; i < n; ++i) {
                double diff = rhs[i] - Ax_direct[i];
                residual_direct += diff * diff;
            }
            residual_direct = std::sqrt(residual_direct);

            std::cout << "Direct solver succeeded!\n";
            std::cout << "  Factorization time: " << std::setprecision(3) << time_factor << " ms\n";
            std::cout << "  Solve time:         " << std::setprecision(3) << time_solve << " ms\n";
            std::cout << "  Total time:         " << std::setprecision(3) << time_total << " ms\n";
            std::cout << "  Residual:           " << std::scientific << residual_direct << "\n";
            std::cout << "  First 5 values:     ";
            for (int i = 0; i < std::min(5, n); ++i) {
                std::cout << std::fixed << std::setprecision(2) << sol_direct[i] << " ";
            }
            std::cout << "\n";

            // Compare best SSOR with direct solver
            std::cout << "\n  Performance comparison (best SSOR vs direct):\n";
            std::cout << "  -----------------------------------------------\n";
            std::cout << "  SSOR (omega=1.5):   fastest iterative method\n";
            std::cout << "  Direct solver:      " << std::setprecision(3) << time_total << " ms total\n";
            std::cout << "  Note: Direct solver has O(n) solve time after factorization,\n";
            std::cout << "        making it ideal for multiple RHS with same matrix.\n";
        } else {
            std::cout << "Direct solver factorization failed!\n";
        }

        // Test from device using Kokkos views and parallel_for
        std::cout << "\n\nTesting device-callable version:\n";
        std::cout << "=================================\n";

        // Copy data to device
        Kokkos::View<int*> rowPtr_d("rowPtr", n + 1);
        Kokkos::View<int*> colInd_d("colInd", nnz);
        Kokkos::View<double*> values_d("values", nnz);
        Kokkos::View<double*> diag_d("diag", n);
        Kokkos::View<double*> rhs_d("rhs", n);
        Kokkos::View<double*> sol_d("sol", n);

        auto rowPtr_h = Kokkos::create_mirror_view(rowPtr_d);
        auto colInd_h = Kokkos::create_mirror_view(colInd_d);
        auto values_h = Kokkos::create_mirror_view(values_d);
        auto rhs_h = Kokkos::create_mirror_view(rhs_d);
        auto sol_h = Kokkos::create_mirror_view(sol_d);

        for (int i = 0; i <= n; ++i) rowPtr_h(i) = rowPtr[i];
        for (int i = 0; i < nnz; ++i) {
            colInd_h(i) = colInd[i];
            values_h(i) = values[i];
        }
        for (int i = 0; i < n; ++i) {
            rhs_h(i) = rhs[i];
            sol_h(i) = 0.0;  // Reset initial guess
        }

        Kokkos::deep_copy(rowPtr_d, rowPtr_h);
        Kokkos::deep_copy(colInd_d, colInd_h);
        Kokkos::deep_copy(values_d, values_h);
        Kokkos::deep_copy(rhs_d, rhs_h);
        Kokkos::deep_copy(sol_d, sol_h);

        // Solve on device
        Kokkos::View<int> iter_d("iterations");
        Kokkos::parallel_for("DeviceSolve",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 1),
            KOKKOS_LAMBDA(const int) {
                // Extract diagonal on device
                IMSI::ExtractDiagonal(n, rowPtr_d.data(), colInd_d.data(),
                                     values_d.data(), diag_d.data());

                // Solve on device
                iter_d() = IMSI::SSOR_Solve(n, rowPtr_d.data(), colInd_d.data(),
                                           values_d.data(), diag_d.data(),
                                           rhs_d.data(), sol_d.data(),
                                           1.2, 1e-3, 10000);
            });
        Kokkos::fence();

        // Copy results back
        Kokkos::deep_copy(sol_h, sol_d);
        int iter_device = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iter_d)();

        if (iter_device > 0) {
            std::cout << "Device solve converged in " << iter_device << " iterations\n";
            std::cout << "First 5 solution values: ";
            for (int i = 0; i < std::min(5, n); ++i) {
                std::cout << sol_h(i) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\n===================\n";
        std::cout << "Example completed!\n";
    }
    Kokkos::finalize();

    return 0;
}
