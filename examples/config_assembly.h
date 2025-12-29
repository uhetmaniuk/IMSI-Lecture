#pragma once

// ============================================================================
// Precision selection: Change 'double' to 'float' to use single precision
// ============================================================================
using Scalar = double;  // Change to 'float' for FP32

// IMSI_COEFF = 1
// \brief Laplace equation for the manufactured solution: u = sin(pi*x) * sin(pi*y)
// IMSI_COEFF = 2
/// \brief Set of functions for 2D oscillatory problem
/// This problem is used in Le Bris, Legoll, Lozinski,
/// "MsFEM a la Crouzeix-Raviart for Highly Oscillatory Elliptic Problems"
/// (2013).
#define IMSI_COEFF 2

#if IMSI_COEFF == 1

// Define material coefficients and forcing term as functors at namespace scope
// These must be KOKKOS_INLINE_FUNCTION compatible for device execution
// NOTE: Must be at namespace scope for CUDA (cannot be local types in main())
struct AxFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    return 1.0;
  }
};

struct AyFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    return 1.0;
  }
};

struct FFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    // Manufactured solution: u = sin(pi*x) * sin(pi*y)
    // => -Laplacian(u) = 2*pi^2 * sin(pi*x) * sin(pi*y)
    constexpr Scalar pi = 3.14159265358979323846;
    using Kokkos::sin;  // Device-safe sin function
    return 2.0 * pi * pi * sin(pi * x) * sin(pi * y);
  }
};

#elif IMSI_COEFF == 2

// Define material coefficients and forcing term as functors at namespace scope
// These must be KOKKOS_INLINE_FUNCTION compatible for device execution
// NOTE: Must be at namespace scope for CUDA (cannot be local types in main())
struct AxFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    using Kokkos::cos;  // Device-safe sin function
    using Kokkos::sin;  // Device-safe sin function
    auto const c = cos(150 * x);
    auto const s = sin(150 * y);
    return Scalar(1 + 100 * c * c * s * s);
  }
};

struct AyFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    using Kokkos::cos;  // Device-safe sin function
    using Kokkos::sin;  // Device-safe sin function
    auto const c = cos(150 * x);
    auto const s = sin(150 * y);
    return Scalar(1 + 100 * c * c * s * s);
  }
};

struct FFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar
  operator()(Scalar x, Scalar y, Scalar z) const
  {
    using Kokkos::sin;  // Device-safe sin function
    return sin(x) * sin(y);
  }
};

#endif
