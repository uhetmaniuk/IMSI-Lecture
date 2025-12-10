#ifndef SSOR_SOLVER_H
#define SSOR_SOLVER_H

#include <Kokkos_Core.hpp>
#include <cmath>

namespace IMSI {

/**
 * @brief Extract diagonal elements from CSR matrix
 *
 * Callable from both host and device.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Output diagonal array (size n, will be filled)
 */
KOKKOS_INLINE_FUNCTION
void ExtractDiagonal(const int n,
                     const int* rowPtr,
                     const int* colInd,
                     const double* values,
                     double* diag)
{
    for (int i = 0; i < n; ++i) {
        diag[i] = 0.0;  // Default if diagonal not found
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            if (colInd[k] == i) {
                diag[i] = values[k];
                break;
            }
        }
    }
}

/**
 * @brief Sequential SSOR solver for sparse linear systems
 *
 * Solves Ax = b where A is a sparse symmetric matrix in CSR format.
 * Uses Symmetric Successive Over-Relaxation iteration.
 *
 * This is a simple sequential implementation callable from both host and device.
 * Designed for small to medium systems (n < 1000) or for use within parallel
 * regions where each thread/team solves an independent system.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Diagonal elements (size n) - call ExtractDiagonal first to fill this
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param omega Relaxation parameter (typically 1.0-1.5, use 1.0 for Gauss-Seidel)
 * @param tol Convergence tolerance for ||b - Ax||_2
 * @param maxIter Maximum number of iterations
 * @return Number of iterations performed (-1 if not converged)
 *
 * Example usage from host:
 * @code
 * int n = 100;
 * int* rowPtr = ...; // CSR format
 * int* colInd = ...;
 * double* values = ...;
 * double* diag = new double[n];
 * double* rhs = ...; // Right-hand side
 * double* sol = ...; // Initial guess
 *
 * IMSI::ExtractDiagonal(n, rowPtr, colInd, values, diag);
 * int iter = IMSI::SSOR_Solve(n, rowPtr, colInd, values, diag,
 *                             rhs, sol, 1.2, 1e-10, 1000);
 * @endcode
 *
 * Example usage from device (TeamPolicy):
 * @code
 * Kokkos::parallel_for("SolveSystems",
 *     Kokkos::TeamPolicy<>(numSystems, Kokkos::AUTO),
 *     KOKKOS_LAMBDA(const member_type& team) {
 *         int isys = team.league_rank();
 *
 *         // Allocate scratch memory for this system
 *         scratch_int_1d rowPtr_s(team.team_scratch(1), n+1);
 *         scratch_int_1d colInd_s(team.team_scratch(1), nnz);
 *         scratch_double_1d values_s(team.team_scratch(1), nnz);
 *         scratch_double_1d diag_s(team.team_scratch(1), n);
 *         scratch_double_1d rhs_s(team.team_scratch(1), n);
 *         scratch_double_1d sol_s(team.team_scratch(1), n);
 *
 *         // ... fill scratch arrays ...
 *
 *         // One thread does the solve (or all threads execute same code)
 *         Kokkos::single(Kokkos::PerTeam(team), [&]() {
 *             IMSI::ExtractDiagonal(n, rowPtr_s.data(), colInd_s.data(),
 *                                   values_s.data(), diag_s.data());
 *             int iter = IMSI::SSOR_Solve(n, rowPtr_s.data(), colInd_s.data(),
 *                                         values_s.data(), diag_s.data(),
 *                                         rhs_s.data(), sol_s.data(),
 *                                         1.2, 1e-10, 100);
 *         });
 *     });
 * @endcode
 */
KOKKOS_INLINE_FUNCTION
int SSOR_Solve(const int n,
               const int* rowPtr,
               const int* colInd,
               const double* values,
               const double* diag,
               const double* rhs,
               double* sol,
               const double omega = 1.0,
               const double tol = 1e-10,
               const int maxIter = 1000)
{
    int iter;
    for (iter = 0; iter < maxIter; ++iter) {
        // Forward sweep (using lower triangular + diagonal)
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                if (colInd[k] != i) {  // Skip diagonal
                    sum += values[k] * sol[colInd[k]];
                }
            }

            if (diag[i] != 0.0) {
                sol[i] = (1.0 - omega) * sol[i] + (omega / diag[i]) * (rhs[i] - sum);
            }
        }

        // Backward sweep (using upper triangular + diagonal)
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                if (colInd[k] != i) {  // Skip diagonal
                    sum += values[k] * sol[colInd[k]];
                }
            }

            if (diag[i] != 0.0) {
                sol[i] = (1.0 - omega) * sol[i] + (omega / diag[i]) * (rhs[i] - sum);
            }
        }

        // Check convergence every 10 iterations
        if ((iter + 1) % 10 == 0) {
            // Compute residual ||b - Ax||_2
            double residual = 0.0;
            for (int i = 0; i < n; ++i) {
                double Ax_i = 0.0;
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                    Ax_i += values[k] * sol[colInd[k]];
                }
                double diff = rhs[i] - Ax_i;
                residual += diff * diff;
            }

            // Use sqrt from cmath
            #ifdef __CUDA_ARCH__
                residual = ::sqrt(residual);
            #else
                residual = std::sqrt(residual);
            #endif

            if (residual < tol) {
                return iter + 1;
            }
        }
    }

    return (iter == maxIter) ? -1 : iter + 1;
}

/**
 * @brief SSOR solver with automatic diagonal extraction
 *
 * Convenience wrapper that extracts diagonal automatically.
 * Note: This allocates a temporary array for diagonal on the stack,
 * so only use for small systems (n < 1000) or ensure sufficient stack space.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param omega Relaxation parameter
 * @param tol Convergence tolerance
 * @param maxIter Maximum iterations
 * @return Number of iterations performed (-1 if not converged)
 */
KOKKOS_INLINE_FUNCTION
int SSOR_Solve_AutoDiag(const int n,
                        const int* rowPtr,
                        const int* colInd,
                        const double* values,
                        const double* rhs,
                        double* sol,
                        const double omega = 1.0,
                        const double tol = 1e-10,
                        const int maxIter = 1000)
{
    // Allocate diagonal on stack (only safe for small n)
    double diag[1024];  // Max system size = 1024
    if (n > 1024) {
        return -1;  // System too large for stack allocation
    }

    ExtractDiagonal(n, rowPtr, colInd, values, diag);
    return SSOR_Solve(n, rowPtr, colInd, values, diag, rhs, sol, omega, tol, maxIter);
}

}  // namespace IMSI

#endif  // SSOR_SOLVER_H
