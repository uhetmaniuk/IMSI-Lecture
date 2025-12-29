#ifndef PCG_SOLVER_H
#define PCG_SOLVER_H

#include <Kokkos_Core.hpp>
#include <cmath>

namespace IMSI {

/**
 * @brief Preconditioned Conjugate Gradient (PCG) solver for sparse linear systems
 *
 * Solves Ax = b where A is a sparse symmetric positive definite matrix in CSR format.
 * Uses diagonal (Jacobi) preconditioning: M = diag(A)
 *
 * This is a simple sequential implementation callable from both host and device.
 * For symmetric positive definite systems, PCG often converges much faster than SSOR.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Diagonal elements (size n) - call ExtractDiagonal first to fill this
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param tol Convergence tolerance for ||r||_2
 * @param maxIter Maximum number of iterations
 * @return Number of iterations performed (-1 if not converged)
 *
 * Algorithm:
 * - r = b - Ax
 * - z = M^{-1} r  (apply preconditioner)
 * - p = z
 * - for k = 0, 1, 2, ...
 *     alpha = (r, z) / (Ap, p)
 *     x += alpha * p
 *     r -= alpha * Ap
 *     if ||r|| < tol: converged
 *     z = M^{-1} r
 *     beta = (r_new, z_new) / (r_old, z_old)
 *     p = z + beta * p
 *
 * Example usage:
 * @code
 * // Prepare matrix and diagonal
 * int n = 100;
 * int* rowPtr, *colInd;
 * double* values, *diag, *rhs, *sol;
 *
 * IMSI::ExtractDiagonal(n, rowPtr, colInd, values, diag);
 *
 * // Solve with PCG
 * int iter = IMSI::PCG_Solve(
 *     n, rowPtr, colInd, values, diag,
 *     rhs, sol, 1e-6, 1000
 * );
 * @endcode
 */
KOKKOS_INLINE_FUNCTION
int
PCG_Solve(
    const int     n,
    const int*    rowPtr,
    const int*    colInd,
    const double* values,
    const double* diag,
    const double* rhs,
    double*       sol,
    const double  tol     = 1e-10,
    const int     maxIter = 1000)
{
  // Allocate workspace on stack (only safe for small n)
  // For larger n, caller should allocate and pass workspace
  const int MAX_STACK_SIZE = 1024;
  if (n > MAX_STACK_SIZE) {
    return -1;  // System too large for stack allocation
  }

  double r[MAX_STACK_SIZE];   // Residual
  double z[MAX_STACK_SIZE];   // Preconditioned residual
  double p[MAX_STACK_SIZE];   // Search direction
  double Ap[MAX_STACK_SIZE];  // A times p

  // Compute initial residual: r = b - Ax
  for (int i = 0; i < n; ++i) {
    double Ax_i = 0.0;
    for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ax_i += values[k] * sol[colInd[k]]; }
    r[i] = rhs[i] - Ax_i;
  }

  // Apply preconditioner: z = M^{-1} r  (diagonal preconditioning)
  for (int i = 0; i < n; ++i) { z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i]; }

  // Initialize search direction: p = z
  for (int i = 0; i < n; ++i) { p[i] = z[i]; }

  // Compute initial (r, z)
  double rz_old = 0.0;
  for (int i = 0; i < n; ++i) { rz_old += r[i] * z[i]; }

  // PCG iteration
  int iter;
  for (iter = 0; iter < maxIter; ++iter) {
    // Compute Ap = A * p
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ap[i] += values[k] * p[colInd[k]]; }
    }

    // Compute (Ap, p)
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) { pAp += Ap[i] * p[i]; }

    // Check for breakdown
    if (pAp <= 0.0) {
      return -1;  // Matrix not positive definite
    }

    // Compute step length: alpha = (r, z) / (Ap, p)
    double alpha = rz_old / pAp;

    // Update solution: x = x + alpha * p
    for (int i = 0; i < n; ++i) { sol[i] += alpha * p[i]; }

    // Update residual: r = r - alpha * Ap
    for (int i = 0; i < n; ++i) { r[i] -= alpha * Ap[i]; }

    // Check convergence every iteration (PCG is cheap per iteration)
    double residual_norm = 0.0;
    for (int i = 0; i < n; ++i) { residual_norm += r[i] * r[i]; }

    if (residual_norm < tol * tol) { return iter + 1; }

    // Apply preconditioner: z = M^{-1} r
    for (int i = 0; i < n; ++i) { z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i]; }

    // Compute new (r, z)
    double rz_new = 0.0;
    for (int i = 0; i < n; ++i) { rz_new += r[i] * z[i]; }

    // Compute beta = (r_new, z_new) / (r_old, z_old)
    double beta = rz_new / rz_old;

    // Update search direction: p = z + beta * p
    for (int i = 0; i < n; ++i) { p[i] = z[i] + beta * p[i]; }

    // Update (r, z) for next iteration
    rz_old = rz_new;
  }

  return (iter == maxIter) ? -1 : iter + 1;
}

/**
 * @brief PCG solver with workspace provided by caller
 *
 * This version doesn't use stack allocation, making it suitable for large systems.
 * Caller must provide workspace arrays.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Diagonal elements (size n)
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param work Workspace array (size 4*n) - will be used internally
 * @param tol Convergence tolerance
 * @param maxIter Maximum iterations
 * @return Number of iterations performed (-1 if not converged)
 */
KOKKOS_INLINE_FUNCTION
int
PCG_Solve_Workspace(
    const int     n,
    const int*    rowPtr,
    const int*    colInd,
    const double* values,
    const double* diag,
    const double* rhs,
    double*       sol,
    double*       work,  // Size 4*n workspace
    const double  tol     = 1e-10,
    const int     maxIter = 1000)
{
  // Partition workspace
  double* r  = work;          // Residual (size n)
  double* z  = work + n;      // Preconditioned residual (size n)
  double* p  = work + 2 * n;  // Search direction (size n)
  double* Ap = work + 3 * n;  // A times p (size n)

  // Compute initial residual: r = b - Ax
  for (int i = 0; i < n; ++i) {
    double Ax_i = 0.0;
    for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ax_i += values[k] * sol[colInd[k]]; }
    r[i] = rhs[i] - Ax_i;
  }

  // Apply preconditioner: z = M^{-1} r
  for (int i = 0; i < n; ++i) { z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i]; }

  // Initialize search direction: p = z
  for (int i = 0; i < n; ++i) { p[i] = z[i]; }

  // Compute initial (r, z)
  double rz_old = 0.0;
  for (int i = 0; i < n; ++i) { rz_old += r[i] * z[i]; }

  // PCG iteration
  int iter;
  for (iter = 0; iter < maxIter; ++iter) {
    // Compute Ap = A * p
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ap[i] += values[k] * p[colInd[k]]; }
    }

    // Compute (Ap, p)
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) { pAp += Ap[i] * p[i]; }

    // Check for breakdown
    if (pAp <= 0.0) {
      return -1;  // Matrix not positive definite
    }

    // Compute step length: alpha = (r, z) / (Ap, p)
    double alpha = rz_old / pAp;

    // Update solution: x = x + alpha * p
    for (int i = 0; i < n; ++i) { sol[i] += alpha * p[i]; }

    // Update residual: r = r - alpha * Ap
    for (int i = 0; i < n; ++i) { r[i] -= alpha * Ap[i]; }

    // Check convergence every 10 iterations
    double residual_norm = 0.0;
    for (int i = 0; i < n; ++i) { residual_norm += r[i] * r[i]; }

    if (residual_norm < tol * tol) { return iter + 1; }

    // Apply preconditioner: z = M^{-1} r
    for (int i = 0; i < n; ++i) { z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i]; }

    // Compute new (r, z)
    double rz_new = 0.0;
    for (int i = 0; i < n; ++i) { rz_new += r[i] * z[i]; }

    // Compute beta = (r_new, z_new) / (r_old, z_old)
    double beta = rz_new / rz_old;

    // Update search direction: p = z + beta * p
    for (int i = 0; i < n; ++i) { p[i] = z[i] + beta * p[i]; }

    // Update (r, z) for next iteration
    rz_old = rz_new;
  }

  return (iter == maxIter) ? -1 : iter + 1;
}

/**
 * @brief Unpreconditioned CG solver (no preconditioning, M = I)
 *
 * Simpler variant without preconditioning. Generally slower convergence
 * than PCG, but useful for well-conditioned problems.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param tol Convergence tolerance
 * @param maxIter Maximum iterations
 * @return Number of iterations performed (-1 if not converged)
 */
KOKKOS_INLINE_FUNCTION
int
CG_Solve(
    const int     n,
    const int*    rowPtr,
    const int*    colInd,
    const double* values,
    const double* rhs,
    double*       sol,
    const double  tol     = 1e-10,
    const int     maxIter = 1000)
{
  // Allocate workspace on stack
  const int MAX_STACK_SIZE = 1024;
  if (n > MAX_STACK_SIZE) { return -1; }

  double r[MAX_STACK_SIZE];   // Residual
  double p[MAX_STACK_SIZE];   // Search direction
  double Ap[MAX_STACK_SIZE];  // A times p

  // Compute initial residual: r = b - Ax
  for (int i = 0; i < n; ++i) {
    double Ax_i = 0.0;
    for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ax_i += values[k] * sol[colInd[k]]; }
    r[i] = rhs[i] - Ax_i;
  }

  // Initialize search direction: p = r
  for (int i = 0; i < n; ++i) { p[i] = r[i]; }

  // Compute initial (r, r)
  double rr_old = 0.0;
  for (int i = 0; i < n; ++i) { rr_old += r[i] * r[i]; }

  // CG iteration
  int iter;
  for (iter = 0; iter < maxIter; ++iter) {
    // Compute Ap = A * p
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ap[i] += values[k] * p[colInd[k]]; }
    }

    // Compute (Ap, p)
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) { pAp += Ap[i] * p[i]; }

    // Check for breakdown
    if (pAp <= 0.0) { return -1; }

    // Compute step length: alpha = (r, r) / (Ap, p)
    double alpha = rr_old / pAp;

    // Update solution: x = x + alpha * p
    for (int i = 0; i < n; ++i) { sol[i] += alpha * p[i]; }

    // Update residual: r = r - alpha * Ap
    for (int i = 0; i < n; ++i) { r[i] -= alpha * Ap[i]; }

    // Check convergence
    double residual_norm = 0.0;
    for (int i = 0; i < n; ++i) { residual_norm += r[i] * r[i]; }

    if (residual_norm < tol * tol) { return iter + 1; }

    // Compute new (r, r)
    double rr_new = 0.0;
    for (int i = 0; i < n; ++i) { rr_new += r[i] * r[i]; }

    // Compute beta = (r_new, r_new) / (r_old, r_old)
    double beta = rr_new / rr_old;

    // Update search direction: p = r + beta * p
    for (int i = 0; i < n; ++i) { p[i] = r[i] + beta * p[i]; }

    // Update (r, r) for next iteration
    rr_old = rr_new;
  }

  return (iter == maxIter) ? -1 : iter + 1;
}

/**
 * @brief PCG solver with SSOR preconditioning
 *
 * Uses SSOR iteration as preconditioner instead of simple diagonal scaling.
 * Applies numSSORSweeps SSOR sweeps to approximately solve M*z = r.
 * Typically 1-2 sweeps are sufficient and more effective than diagonal preconditioning.
 *
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Diagonal elements (size n)
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param work Workspace array (size 4*n)
 * @param omega SSOR relaxation parameter (typically 1.0-1.5)
 * @param numSSORSweeps Number of SSOR sweeps per preconditioner application (typically 1-2)
 * @param tol Convergence tolerance
 * @param maxIter Maximum PCG iterations
 * @return Number of iterations performed (-1 if not converged)
 *
 * Performance note:
 * - More expensive per iteration than diagonal preconditioning
 * - But may reduce total iterations significantly for ill-conditioned systems
 * - Cost per iteration: ~numSSORSweeps * (2 SpMV) + regular PCG cost
 *
 * Example usage:
 * @code
 * double work[4*n];
 * int iter = IMSI::PCG_Solve_SSOR_Precond(
 *     n, rowPtr, colInd, values, diag,
 *     rhs, sol, work,
 *     1.2,    // omega for SSOR
 *     1,      // 1 SSOR sweep per preconditioner application
 *     1e-6, 1000
 * );
 * @endcode
 */
KOKKOS_INLINE_FUNCTION
int
PCG_Solve_SSOR_Precond(
    const int     n,
    const int*    rowPtr,
    const int*    colInd,
    const double* values,
    const double* diag,
    const double* rhs,
    double*       sol,
    double*       work,  // Size 4*n workspace
    const double  omega         = 1.0,
    const int     numSSORSweeps = 1,
    const double  tol           = 1e-10,
    const int     maxIter       = 1000)
{
  // Partition workspace
  double* r  = work;          // Residual (size n)
  double* z  = work + n;      // Preconditioned residual (size n)
  double* p  = work + 2 * n;  // Search direction (size n)
  double* Ap = work + 3 * n;  // A times p (size n)

  // Compute initial residual: r = b - Ax
  for (int i = 0; i < n; ++i) {
    double Ax_i = 0.0;
    for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ax_i += values[k] * sol[colInd[k]]; }
    r[i] = rhs[i] - Ax_i;
  }

  // Apply SSOR preconditioner: z = M^{-1} r
  // Solve M*z = r approximately using SSOR iteration
  // Start with z = 0
  for (int i = 0; i < n; ++i) { z[i] = 0.0; }

  // Apply numSSORSweeps SSOR sweeps
  for (int sweep = 0; sweep < numSSORSweeps; ++sweep) {
    // Forward sweep
    for (int i = 0; i < n; ++i) {
      double sum = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
        if (colInd[k] != i) { sum += values[k] * z[colInd[k]]; }
      }
      if (diag[i] != 0.0) { z[i] = (1.0 - omega) * z[i] + (omega / diag[i]) * (r[i] - sum); }
    }

    // Backward sweep
    for (int i = n - 1; i >= 0; --i) {
      double sum = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
        if (colInd[k] != i) { sum += values[k] * z[colInd[k]]; }
      }
      if (diag[i] != 0.0) { z[i] = (1.0 - omega) * z[i] + (omega / diag[i]) * (r[i] - sum); }
    }
  }

  // Initialize search direction: p = z
  for (int i = 0; i < n; ++i) { p[i] = z[i]; }

  // Compute initial (r, z)
  double rz_old = 0.0;
  for (int i = 0; i < n; ++i) { rz_old += r[i] * z[i]; }

  // PCG iteration
  int iter;
  for (iter = 0; iter < maxIter; ++iter) {
    // Compute Ap = A * p
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ap[i] += values[k] * p[colInd[k]]; }
    }

    // Compute (Ap, p)
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) { pAp += Ap[i] * p[i]; }

    // Check for breakdown
    if (pAp <= 0.0) {
      return -1;  // Matrix not positive definite
    }

    // Compute step length: alpha = (r, z) / (Ap, p)
    double alpha = rz_old / pAp;

    // Update solution: x = x + alpha * p
    for (int i = 0; i < n; ++i) { sol[i] += alpha * p[i]; }

    // Update residual: r = r - alpha * Ap
    for (int i = 0; i < n; ++i) { r[i] -= alpha * Ap[i]; }

    // Check convergence
    double residual_norm = 0.0;
    for (int i = 0; i < n; ++i) { residual_norm += r[i] * r[i]; }

    if (residual_norm < tol * tol) { return iter + 1; }

    // Apply SSOR preconditioner: z = M^{-1} r
    // Reset z to 0 and apply SSOR sweeps
    for (int i = 0; i < n; ++i) { z[i] = 0.0; }

    for (int sweep = 0; sweep < numSSORSweeps; ++sweep) {
      // Forward sweep
      for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
          if (colInd[k] != i) { sum += values[k] * z[colInd[k]]; }
        }
        if (diag[i] != 0.0) { z[i] = (1.0 - omega) * z[i] + (omega / diag[i]) * (r[i] - sum); }
      }

      // Backward sweep
      for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
          if (colInd[k] != i) { sum += values[k] * z[colInd[k]]; }
        }
        if (diag[i] != 0.0) { z[i] = (1.0 - omega) * z[i] + (omega / diag[i]) * (r[i] - sum); }
      }
    }

    // Compute new (r, z)
    double rz_new = 0.0;
    for (int i = 0; i < n; ++i) { rz_new += r[i] * z[i]; }

    // Compute beta = (r_new, z_new) / (r_old, z_old)
    double beta = rz_new / rz_old;

    // Update search direction: p = z + beta * p
    for (int i = 0; i < n; ++i) { p[i] = z[i] + beta * p[i]; }

    // Update (r, z) for next iteration
    rz_old = rz_new;
  }

  return (iter == maxIter) ? -1 : iter + 1;
}

/**
 * @brief Team-parallel PCG solver with Jacobi preconditioning for Kokkos teams
 *
 * This version parallelizes all PCG operations within a Kokkos team using TeamThreadRange.
 * Designed for GPU execution where each team solves one linear system.
 *
 * All team members participate in:
 * - SpMV (sparse matrix-vector product)
 * - Dot products (parallel reductions)
 * - Vector updates (AXPY operations)
 *
 * @tparam TeamMember Kokkos team member type
 * @param team The team member handle
 * @param n Matrix dimension
 * @param rowPtr CSR row pointers (size n+1)
 * @param colInd CSR column indices
 * @param values CSR values
 * @param diag Diagonal elements for Jacobi preconditioning (size n)
 * @param rhs Right-hand side (size n)
 * @param sol Solution vector (size n) - input: initial guess, output: solution
 * @param work Workspace array (size 4*n): partitioned as [r, z, p, Ap]
 * @param tol Convergence tolerance for ||r||_2
 * @param maxIter Maximum number of iterations
 * @return Number of iterations performed (-1 if not converged)
 */
template <typename TeamMember>
KOKKOS_INLINE_FUNCTION int
PCG_Solve_Team(
    const TeamMember& team,
    const int         n,
    const int*        rowPtr,
    const int*        colInd,
    const double*     values,
    const double*     diag,
    const double*     rhs,
    double*           sol,
    double*           work,
    const double      tol     = 1e-10,
    const int         maxIter = 1000)
{
  // Partition workspace
  double* r  = work;          // Residual (size n)
  double* z  = work + n;      // Preconditioned residual (size n)
  double* p  = work + 2 * n;  // Search direction (size n)
  double* Ap = work + 3 * n;  // A times p (size n)

  // Compute initial residual: r = b - Ax (team-parallel SpMV)
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) {
    double Ax_i = 0.0;
    for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { Ax_i += values[k] * sol[colInd[k]]; }
    r[i] = rhs[i] - Ax_i;
  });
  team.team_barrier();

  // Apply Jacobi preconditioner: z = M^{-1} r
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) {
    z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i];
  });
  team.team_barrier();

  // Initialize search direction: p = z
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) { p[i] = z[i]; });
  team.team_barrier();

  // Compute initial (r, z) using team-parallel reduction
  double rz_old = 0.0;
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team, n), [&](int i, double& sum) { sum += r[i] * z[i]; }, rz_old);
  team.team_barrier();

  // PCG iteration
  int iter;
  for (iter = 0; iter < maxIter; ++iter) {
    // Compute Ap = A * p (team-parallel SpMV)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) {
      double sum = 0.0;
      for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) { sum += values[k] * p[colInd[k]]; }
      Ap[i] = sum;
    });
    team.team_barrier();

    // Compute (Ap, p) using team-parallel reduction
    double pAp = 0.0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, n), [&](int i, double& sum) { sum += Ap[i] * p[i]; }, pAp);
    team.team_barrier();

    // Check for breakdown
    if (pAp <= 0.0) { return -1; }

    // Compute step length: alpha = (r, z) / (Ap, p)
    double alpha = rz_old / pAp;

    // Update solution and residual in parallel
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) {
      sol[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    });
    team.team_barrier();

    // Check convergence using team-parallel reduction
    double residual_norm_sq = 0.0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, n),
        [&](int i, double& sum) { sum += r[i] * r[i]; },
        residual_norm_sq);
    team.team_barrier();

    if (residual_norm_sq < tol * tol) { return iter + 1; }

    // Apply Jacobi preconditioner: z = M^{-1} r
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) {
      z[i] = (diag[i] != 0.0) ? r[i] / diag[i] : r[i];
    });
    team.team_barrier();

    // Compute new (r, z) using team-parallel reduction
    double rz_new = 0.0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, n), [&](int i, double& sum) { sum += r[i] * z[i]; }, rz_new);
    team.team_barrier();

    // Compute beta = (r_new, z_new) / (r_old, z_old)
    double beta = rz_new / rz_old;

    // Update search direction: p = z + beta * p
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [&](int i) { p[i] = z[i] + beta * p[i]; });
    team.team_barrier();

    // Update (r, z) for next iteration
    rz_old = rz_new;
  }

  return (iter == maxIter) ? -1 : iter + 1;
}

}  // namespace IMSI

#endif  // PCG_SOLVER_H
