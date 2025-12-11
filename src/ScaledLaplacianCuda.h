#pragma once

#include <functional>
#include <optional>

#include "../main_config.h"
#include "Element.h"
#include "MathUtils.h"  // For InverseInPlace
#include "MathUtilsCuda.h"
#include "MeshUtils.h"
#include "QuadratureRule.h"
#include "SparseMatrix.hpp"
#include "SymmetricSparse.hpp"
#include "fe2DQ1.h"  // For host-side MFEM assembly
#include "fe2DQ1Cuda.h"

namespace IMSI {

/// \brief Device-compatible binary search for sorted array
/// Returns index of found element or position where it would be inserted
template <typename T>
KOKKOS_INLINE_FUNCTION int
lower_bound_device(const T* array, int size, T value)
{
  int left  = 0;
  int right = size;
  while (left < right) {
    int mid = left + (right - left) / 2;
    if (array[mid] < value) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  return left;
}

/// \brief Functor for MFEM element assembly
/// Must be at namespace scope for CUDA compatibility
/// \tparam ExecutionSpace Kokkos execution space (e.g., Kokkos::Cuda, Kokkos::OpenMP)
/// \tparam FuncX Functor type for ax coefficient
/// \tparam FuncY Functor type for ay coefficient
/// \tparam FuncF Functor type for f coefficient
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
struct MFEMAssemblyFunctor
{
  Kokkos::View<int*, ExecutionSpace>    eleList_device;
  Kokkos::View<int*, ExecutionSpace>    cellTypes;
  Kokkos::View<double*, ExecutionSpace> nodeCoords;
  Kokkos::View<int**, ExecutionSpace>   cellToNode;
  int                                   numEle;

  // Quadrature data
  Kokkos::View<double*, ExecutionSpace> quadWeight;
  Kokkos::View<double*, ExecutionSpace> quadXi;
  Kokkos::View<double*, ExecutionSpace> quadEta;
  Kokkos::View<double*, ExecutionSpace> quadZeta;
  int                                   ruleLen;

  // Coefficient functors
  FuncX ax_func;
  FuncY ay_func;
  FuncF f_func;

  static int const ratio            = 32;
  static int constexpr numFineNodes = (ratio + 1) * (ratio + 1);
  static int constexpr numFineEle   = ratio * ratio;

  typedef typename Kokkos::TeamPolicy<ExecutionSpace>::member_type team_member;
  typedef typename ExecutionSpace::scratch_memory_space            scratch_space;
  typedef Kokkos::View<int*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>    scratch_int_1d;
  typedef Kokkos::View<double*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_double_1d;

  KOKKOS_FUNCTION
  void
  operator()(const team_member& teamMember) const
  {
    // Get element index
    const int ieleCoarse = teamMember.league_rank();

    // Allocate scratch pad memory at team level (level 1)
    // 2 vectors of size numFineNodes
    scratch_double_1d rhs(teamMember.team_scratch(1), numFineNodes);

    // Extract globalToFree and freeToGlobal mappings for interior DOFs
    // Interior DOFs are those where (ix > 0) && (ix < ratio) && (iy > 0) && (iy < ratio)
    // Number of free DOFs = (ratio-1) * (ratio-1)
    int const           numFreeDofs = (ratio - 1) * (ratio - 1);
    scratch_int_1d      globalToFree(teamMember.team_scratch(1), numFineNodes);
    scratch_int_1d      freeToGlobal(teamMember.team_scratch(1), numFreeDofs);

    // Extract globalToBoundary and boundaryToGlobal mappings for boundary DOFs
    // Boundary DOFs are those where (ix == 0) || (ix == ratio) || (iy == 0) || (iy == ratio)
    int const      numBoundaryDofs = numFineNodes - numFreeDofs;  // 128 for ratio=32
    scratch_int_1d globalToBoundary(teamMember.team_scratch(1), numFineNodes);
    scratch_int_1d boundaryToGlobal(teamMember.team_scratch(1), numBoundaryDofs);

    // Setup basis functions on boundaries
    // phi is a matrix of size (numVectors x numFineNodes) stored in column-major order
    // Each column represents one basis function evaluated at all fine nodes
    scratch_double_1d phi(teamMember.team_scratch(1), numVectors * numFineNodes);

    // Initialize phi to zero
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numVectors * numFineNodes), [&](int i) {
      phi(i) = 0.0;
    });

    // Initialize rhs to 0, globalToFree to -1 (indicating constrained/boundary DOF)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFineNodes), [&](int i) {
      rhs(i)      = 0;
      globalToFree(i) = -1;
      globalToBoundary(i) = -1;
    });

    // 8 vectors for static condensation: 4 RHS + 4 solutions (one per basis function)
    constexpr int     numVectors = 4;  // Number of coarse element nodes
    scratch_double_1d btmp(teamMember.team_scratch(1), numFreeDofs * numVectors);  // RHS vectors
    scratch_double_1d utmp(teamMember.team_scratch(1), numFreeDofs * numVectors);  // Solution vectors
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs * numVectors), [&](int i) {
      btmp(i) = 0.0;
      utmp(i) = 0.0;
    });

    // Set a barrier as we will update phi, globalToFree, freeToGlobal, boundaryToGlobal, globalToBoundary
    teamMember.team_barrier();

    // Build DOF mappings and corner basis functions (serial operations)
    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
      // Build both interior and boundary mappings
      int freeCount     = 0;
      int boundaryCount = 0;
      for (int iy = 0; iy <= ratio; ++iy) {
        for (int ix = 0; ix <= ratio; ++ix) {
          int const nodeID = ix + iy * (ratio + 1);
          // Check if interior or boundary
          if ((ix > 0) && (ix < ratio) && (iy > 0) && (iy < ratio)) {
            // Interior node
            globalToFree(nodeID)       = freeCount;
            freeToGlobal(freeCount)    = nodeID;
            freeCount += 1;
          } else {
            // Boundary node
            globalToBoundary(nodeID)     = boundaryCount;
            boundaryToGlobal(boundaryCount) = nodeID;
            boundaryCount += 1;
          }
        }
      }
      // Handle 4 corners for phi basis functions (avoids race conditions)
      // Corner (0, 0) - Basis 0
      phi(0 + 0 * numFineNodes) = 1.0;
      // Corner (ratio, 0) - Basis 1
      int in = ratio;
      phi(in + 1 * numFineNodes) = 1.0;
      // Corner (0, ratio) - Basis 3
      in = ratio * (ratio + 1);
      phi(in + 3 * numFineNodes) = 1.0;
      // Corner (ratio, ratio) - Basis 2
      in = ratio + ratio * (ratio + 1);
      phi(in + 2 * numFineNodes) = 1.0;
    });

    // Build boundary basis functions in parallel (skip corners)
    // Left and right edges (ix = 0 and ix = ratio)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 1, ratio), [&](int is) {
      double const s = double(is) / double(ratio);  // abscissa along an edge
      // Left edge (ix = 0)
      int in = is * (ratio + 1);
      phi(in + 0 * numFineNodes) = 1.0 - s;  // Basis function 0
      phi(in + 3 * numFineNodes) = s;        // Basis function 3
      // Right edge (ix = ratio)
      in = ratio + is * (ratio + 1);
      phi(in + 1 * numFineNodes) = 1.0 - s;  // Basis function 1
      phi(in + 2 * numFineNodes) = s;        // Basis function 2
    });

    // Bottom and top edges (iy = 0 and iy = ratio, skip corners)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 1, ratio), [&](int is) {
      double const s = double(is) / double(ratio);  // abscissa along an edge
      // Bottom edge (iy = 0)
      int in = is;
      phi(in + 0 * numFineNodes) = 1.0 - s;  // Basis function 0
      phi(in + 1 * numFineNodes) = s;        // Basis function 1
      // Top edge (iy = ratio)
      in = is + ratio * (ratio + 1);
      phi(in + 3 * numFineNodes) = 1.0 - s;  // Basis function 3
      phi(in + 2 * numFineNodes) = s;        // Basis function 2
    });

    // Build K_ii: CSR matrix for interior DOFs only
    // CSR matrix row pointer (we'll compute actual nnz before allocating colIdx and values)
    scratch_int_1d matRowPtr_ii(teamMember.team_scratch(1), numFreeDofs + 1);

    // Get the sparsity graph of K_ii (interior-interior stiffness matrix)
    // Map from free DOF index to sparsity pattern
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int iFree, int& update, bool final) {
      int const iGlobal = freeToGlobal(iFree);
      int const ix      = iGlobal % (ratio + 1);
      int const iy      = iGlobal / (ratio + 1);

      // Count interior neighbors only
      int count = 1;  // Diagonal
      int const hasWest  = (ix > 1);
      int const hasEast  = (ix < ratio - 1);
      int const hasSouth = (iy > 1);
      int const hasNorth = (iy < ratio - 1);

      // South neighbors (iy - 1 row)
      if (hasSouth) {
        count += hasWest;   // SW
        count++;            // S (always present when hasSouth)
        count += hasEast;   // SE
      }
      // Center row neighbors
      count += hasWest + hasEast;  // W and E
      // North neighbors (iy + 1 row)
      if (hasNorth) {
        count += hasWest;   // NW
        count++;            // N (always present when hasNorth)
        count += hasEast;   // NE
      }

      if (final) {
        matRowPtr_ii(iFree) = update;
        if (iFree == numFreeDofs - 1) matRowPtr_ii(numFreeDofs) = update + count;
      }
      update += count;
    });
    teamMember.team_barrier();

    // Now allocate colIdx and values for K_ii with the exact number of non-zeros
    int const         actualNnz_ii = matRowPtr_ii(numFreeDofs);
    scratch_int_1d    matColIdx_ii(teamMember.team_scratch(1), actualNnz_ii);
    scratch_double_1d matValues_ii(teamMember.team_scratch(1), actualNnz_ii);

    // Build column indices for K_ii sparsity pattern (in free DOF numbering)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int iFree) {
      int const iGlobal = freeToGlobal(iFree);
      int const ix      = iGlobal % (ratio + 1);
      int const iy      = iGlobal / (ratio + 1);
      int       offset  = matRowPtr_ii(iFree);

      // Add column indices for interior neighbors (in free DOF numbering)
      // South neighbors (iy - 1 row)
      if (iy > 1) {
        // SW neighbor
        if (ix > 1) {
          int jGlobal = iGlobal - 1 - (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
        // S neighbor
        {
          int jGlobal = iGlobal - (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
        // SE neighbor
        if (ix < ratio - 1) {
          int jGlobal = iGlobal + 1 - (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
      }
      // W neighbor
      if (ix > 1) {
        int jGlobal = iGlobal - 1;
        int jFree   = globalToFree(jGlobal);
        if (jFree != -1) matColIdx_ii(offset++) = jFree;
      }
      // Diagonal
      matColIdx_ii(offset++) = iFree;
      // E neighbor
      if (ix < ratio - 1) {
        int jGlobal = iGlobal + 1;
        int jFree   = globalToFree(jGlobal);
        if (jFree != -1) matColIdx_ii(offset++) = jFree;
      }
      // North neighbors (iy + 1 row)
      if (iy < ratio - 1) {
        // NW neighbor
        if (ix > 1) {
          int jGlobal = iGlobal - 1 + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
        // N neighbor
        {
          int jGlobal = iGlobal + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
        // NE neighbor
        if (ix < ratio - 1) {
          int jGlobal = iGlobal + 1 + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1) matColIdx_ii(offset++) = jFree;
        }
      }
    });

    // Initialize K_ii matrix values to zero
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, actualNnz_ii), [&](int i) {
      matValues_ii(i) = 0.0;
    });
    teamMember.team_barrier();

    // Build K_b: CSR matrix for boundary DOFs (numBoundaryDofs x numFineNodes)
    // K_b = [K_bi, K_bb] includes both boundary-interior and boundary-boundary couplings
    // We store K_b in global node numbering for columns (easier indexing)
    scratch_int_1d matRowPtr_b(teamMember.team_scratch(1), numBoundaryDofs + 1);

    // Get the sparsity graph of K_b (boundary rows, all columns)
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(teamMember, numBoundaryDofs), [&](int iBoundary, int& update, bool final) {
      int const iGlobal = boundaryToGlobal(iBoundary);
      int const ix      = iGlobal % (ratio + 1);
      int const iy      = iGlobal / (ratio + 1);

      // Count all neighbors (interior and boundary) in global numbering
      int count = 1;  // Diagonal
      int const hasWest = (ix > 0);
      int const hasEast = (ix < ratio);
      count += hasWest + hasEast;                          // W and E neighbors
      if (iy > 0) { count += 1 + hasWest + hasEast; }     // S, SW, SE neighbors
      if (iy < ratio) { count += 1 + hasWest + hasEast; } // N, NW, NE neighbors

      if (final) {
        matRowPtr_b(iBoundary) = update;
        if (iBoundary == numBoundaryDofs - 1) matRowPtr_b(numBoundaryDofs) = update + count;
      }
      update += count;
    });

    teamMember.team_barrier();

    // Now allocate colIdx and values for K_b with the exact number of non-zeros
    int const         actualNnz_b = matRowPtr_b(numBoundaryDofs);
    scratch_int_1d    matColIdx_b(teamMember.team_scratch(1), actualNnz_b);
    scratch_double_1d matValues_b(teamMember.team_scratch(1), actualNnz_b);

    // Build column indices for K_b sparsity pattern (in global node numbering)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numBoundaryDofs), [&](int iBoundary) {
      int const iGlobal = boundaryToGlobal(iBoundary);
      int const ix      = iGlobal % (ratio + 1);
      int const iy      = iGlobal / (ratio + 1);
      int       offset  = matRowPtr_b(iBoundary);

      // Add column indices for all neighbors (in global node numbering)
      if (iy > 0) {
        if (ix > 0) matColIdx_b(offset++) = iGlobal - 1 - (ratio + 1);
        matColIdx_b(offset++) = iGlobal - (ratio + 1);
        if (ix < ratio) matColIdx_b(offset++) = iGlobal + 1 - (ratio + 1);
      }
      if (ix > 0) matColIdx_b(offset++) = iGlobal - 1;
      matColIdx_b(offset++) = iGlobal;  // Diagonal
      if (ix < ratio) matColIdx_b(offset++) = iGlobal + 1;
      if (iy < ratio) {
        if (ix > 0) matColIdx_b(offset++) = iGlobal - 1 + (ratio + 1);
        matColIdx_b(offset++) = iGlobal + (ratio + 1);
        if (ix < ratio) matColIdx_b(offset++) = iGlobal + 1 + (ratio + 1);
      }
    });
    teamMember.team_barrier();

    // Initialize K_b matrix values to zero
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, actualNnz_b), [&](int i) {
      matValues_b(i) = 0.0;
    });
    teamMember.team_barrier();

    // Get coarse element coordinates from global mesh
    int const     eleID = eleList_device(ieleCoarse);
    constexpr int nNodesCoarse = 4;  // Q1 coarse element
    constexpr int sdim         = 2;

    // Extract coarse element node coordinates
    double coords_coarse[nNodesCoarse * sdim];
    for (int i = 0; i < nNodesCoarse; ++i) {
      int const nodeID           = cellToNode(eleID, i);
      coords_coarse[i * sdim + 0] = nodeCoords(nodeID * sdim + 0);
      coords_coarse[i * sdim + 1] = nodeCoords(nodeID * sdim + 1);
    }

    // Compute fine mesh spacing (assuming Q1 coarse element: 4 nodes)
    double const hx = (coords_coarse[2] - coords_coarse[0]) / double(ratio);
    double const hy = (coords_coarse[7] - coords_coarse[1]) / double(ratio);

    // Assemble fine elements (scalar path only, no SIMD)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFineEle), [&](int ieleLocal) {
      int const ix = ieleLocal % ratio;
      int const iy = ieleLocal / ratio;

      // Node list for this fine element (local Q1 element)
      constexpr int nNodesFine = 4;
      int nodeList[nNodesFine];
      nodeList[0] = ix + iy * (ratio + 1);
      nodeList[1] = ix + 1 + iy * (ratio + 1);
      nodeList[2] = ix + 1 + (iy + 1) * (ratio + 1);
      nodeList[3] = ix + (iy + 1) * (ratio + 1);

      // Compute fine element coordinates
      double coords_fine[nNodesFine * sdim];
      coords_fine[0] = coords_coarse[0] + ix * hx;
      coords_fine[1] = coords_coarse[1] + iy * hy;
      coords_fine[2] = coords_fine[0] + hx;
      coords_fine[3] = coords_fine[1];
      coords_fine[4] = coords_fine[2];
      coords_fine[5] = coords_fine[3] + hy;
      coords_fine[6] = coords_fine[0];
      coords_fine[7] = coords_fine[5];

      // Allocate element arrays
      double rFineEle[nNodesFine];
      double kFineEle[nNodesFine * nNodesFine];

      // Call element assembly (note: need to use template syntax for static member function)
      ScaledLaplacianCuda<ExecutionSpace, FuncX, FuncY, FuncF>::template ElementaryDataQ1<double>(
          coords_fine,
          &quadWeight(0),
          &quadXi(0),
          &quadEta(0),
          &quadZeta(0),
          ruleLen,
          ax_func,
          ay_func,
          f_func,
          rFineEle,
          kFineEle);

      // Scatter RHS to local system (atomic adds for race-free assembly)
      for (int in = 0; in < nNodesFine; ++in) {
        Kokkos::atomic_add(&rhs(nodeList[in]), rFineEle[in]);
      }

      // Scatter stiffness matrix to K_ii, K_b, and btmp in a single pass
      for (int in = 0; in < nNodesFine; ++in) {
        int const iGlobal = nodeList[in];
        int const iFree   = globalToFree(iGlobal);

        for (int jn = 0; jn < nNodesFine; ++jn) {
          int const    jGlobal = nodeList[jn];
          double const k_val   = kFineEle[in + jn * nNodesFine];

          if (iFree != -1) {
            // Interior row
            int const jFree = globalToFree(jGlobal);
            if (jFree != -1) {
              // Interior-interior coupling: add to K_ii
              int const rowBegin = matRowPtr_ii(iFree);
              int const rowEnd   = matRowPtr_ii(iFree + 1);
              int const pos      = lower_bound_device(&matColIdx_ii(rowBegin), rowEnd - rowBegin, jFree);
              Kokkos::atomic_add(&matValues_ii(rowBegin + pos), k_val);
            } else {
              // Interior-boundary coupling: add to btmp (RHS for static condensation)
              for (int iv = 0; iv < numVectors; ++iv) {
                double const phi_val = phi(jGlobal + iv * numFineNodes);
                Kokkos::atomic_add(&btmp(iFree + iv * numFreeDofs), -k_val * phi_val);
              }
            }
          } else {
            // Boundary row: add to K_b (both K_bi and K_bb)
            int const iBoundary = globalToBoundary(iGlobal);
            int const rowBegin  = matRowPtr_b(iBoundary);
            int const rowEnd    = matRowPtr_b(iBoundary + 1);
            int const pos       = lower_bound_device(&matColIdx_b(rowBegin), rowEnd - rowBegin, jGlobal);
            Kokkos::atomic_add(&matValues_b(rowBegin + pos), k_val);
          }
        }
      }
    });

    teamMember.team_barrier();

    //  Solve the linear system
    //  Compute Phi^T K Phi
    //  Compute Phi^T rhs
    //  Scatter them into the coarse matrix and right hand side
  }
};

/// \brief CUDA ScaledLaplacian class for 2D problems
///
/// This class provides CUDA-accelerated assembly for the scaled Laplacian operator:
///   -∇·(α∇u) = f
///
/// Supported elements:
///   - Q1 (bilinear quadrilateral) in 2D
///   - MFEM_L (Multiscale FEM with static condensation) in 2D
///
/// Note: Q2, 1D, and 3D cases are not implemented (per user requirements)
///
/// Template parameters:
///   - ExecutionSpace: Kokkos execution space (e.g., Kokkos::Cuda, Kokkos::OpenMP)
///   - FuncX, FuncY, FuncF: Functor types for coefficients ax, ay, and f
///     Each must have operator()(double x, double y, double z) marked KOKKOS_INLINE_FUNCTION
///
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
class ScaledLaplacianCuda
{
 public:
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  ScaledLaplacianCuda(
      const MeshConnectivity<>& meshData,
      RuleType                  quadRule,
      int                       quadOrder,
      FuncX                     ax_in,
      FuncY                     ay_in,
      FuncF                     f_in)
      : meshInfo(meshData), ruleType(quadRule), ruleOrder(quadOrder), ax_func(ax_in), ay_func(ay_in), f_func(f_in)
  {
    auto const sdim = meshInfo.mesh.GetSpatialDimension();
    if (sdim != 2) { throw std::runtime_error("ScaledLaplacianCuda is only implemented for 2D problems"); }

    getQuadrature(ruleType, sdim, ruleOrder, ruleLength, weight, xi, eta, zeta);

    // Copy quadrature data to device
    quadWeight_d = Kokkos::View<double*, ExecutionSpace>("quadWeight", ruleLength);
    quadXi_d     = Kokkos::View<double*, ExecutionSpace>("quadXi", ruleLength);
    quadEta_d    = Kokkos::View<double*, ExecutionSpace>("quadEta", ruleLength);
    quadZeta_d   = Kokkos::View<double*, ExecutionSpace>("quadZeta", ruleLength);

    auto quadWeight_h = Kokkos::create_mirror_view(quadWeight_d);
    auto quadXi_h     = Kokkos::create_mirror_view(quadXi_d);
    auto quadEta_h    = Kokkos::create_mirror_view(quadEta_d);
    auto quadZeta_h   = Kokkos::create_mirror_view(quadZeta_d);

    for (int i = 0; i < ruleLength; ++i) {
      quadWeight_h(i) = weight[i];
      quadXi_h(i)     = xi[i];
      quadEta_h(i)    = eta[i];
      quadZeta_h(i)   = zeta[i];
    }

    Kokkos::deep_copy(quadWeight_d, quadWeight_h);
    Kokkos::deep_copy(quadXi_d, quadXi_h);
    Kokkos::deep_copy(quadEta_d, quadEta_h);
    Kokkos::deep_copy(quadZeta_d, quadZeta_h);
  }

  /// \brief Assemble linear system on CUDA using graph coloring
  ///
  /// \param[out] rhs Right-hand side vector (device memory)
  /// \param[in] matRowPtr CSR row pointer (device memory)
  /// \param[in] matColIdx CSR column indices (device memory)
  /// \param[out] matValues CSR values (device memory)
  ///
  /// Uses graph coloring to enable conflict-free parallel assembly on GPU
  void
  GetLinearSystem(
      Kokkos::View<double*, ExecutionSpace> rhs,
      Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
      Kokkos::View<int*, ExecutionSpace>    matColIdx,
      Kokkos::View<double*, ExecutionSpace> matValues);

  void
  GetLinearSystemMFEM(
      Kokkos::View<double*, ExecutionSpace> rhs,
      Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
      Kokkos::View<int*, ExecutionSpace>    matColIdx,
      Kokkos::View<double*, ExecutionSpace> matValues);

  int const ratio = 32;  // For MFEM_L fine mesh refinement

  /// \brief Process MFEM_L elements on host (requires sparse solve)
  ///
  /// MFEM_L elements require local sparse solves, which are done on CPU
  /// using the SymmetricSparse solver. This function handles MFEM_L elements
  /// separately from Q1 elements.
  ///
  void
  ProcessMFEMElements(
      Kokkos::View<double*, HostSpace> rhs_h,
      Kokkos::View<size_t*, HostSpace> matRowPtr_h,
      Kokkos::View<int*, HostSpace>    matColIdx_h,
      Kokkos::View<double*, HostSpace> matValues_h);

 protected:
  const MeshConnectivity<> meshInfo;

  /// Coefficient functors (stored as member variables)
  FuncX ax_func;
  FuncY ay_func;
  FuncF f_func;

  RuleType ruleType   = RuleType::Gauss;
  int      ruleOrder  = 1;
  int      ruleLength = 0;

  std::vector<double> weight;
  std::vector<double> xi, eta, zeta;

  // Device copies of quadrature data
  Kokkos::View<double*, ExecutionSpace> quadWeight_d;
  Kokkos::View<double*, ExecutionSpace> quadXi_d;
  Kokkos::View<double*, ExecutionSpace> quadEta_d;
  Kokkos::View<double*, ExecutionSpace> quadZeta_d;

  /// \brief Structure to hold mesh data on device
  struct MeshDeviceData
  {
    Kokkos::View<int*, ExecutionSpace>    cellTypes;
    Kokkos::View<double*, ExecutionSpace> nodeCoords;
    Kokkos::View<int**, ExecutionSpace>   cellToNode;
  };

  /// \brief Copy mesh data to device
  ///
  /// Extracts mesh topology and geometry to device-accessible Kokkos::Views
  /// \return MeshDeviceData structure with device views
  MeshDeviceData
  CopyMeshToDevice() const
  {
    int const numCells = meshInfo.mesh.NumberCells();
    int const numNodes = meshInfo.mesh.NumberVertices();
    int const sdim     = meshInfo.mesh.GetSpatialDimension();

    MeshDeviceData data;

    // Copy cell types
    data.cellTypes = Kokkos::View<int*, ExecutionSpace>("cellTypes", numCells);
    {
      auto  cellTypes_h   = Kokkos::create_mirror_view(data.cellTypes);
      auto* cellTypes_ptr = meshInfo.mesh.GetCellType().data();
      Kokkos::parallel_for("CellTypes_Copy", Kokkos::RangePolicy<HostSpace>(0, numCells), [=](const int ic) {
        cellTypes_h(ic) = static_cast<int>(cellTypes_ptr[ic]);
      });
      Kokkos::deep_copy(data.cellTypes, cellTypes_h);
    }

    // Copy node coordinates (interleaved)
    data.nodeCoords = Kokkos::View<double*, ExecutionSpace>("nodeCoords", numNodes * sdim);
    {
      auto  nodeCoords_h = Kokkos::create_mirror_view(data.nodeCoords);
      auto& mesh_ref     = meshInfo.mesh;
      Kokkos::parallel_for("NodeCoords_Copy", Kokkos::RangePolicy<HostSpace>(0, numNodes), [=, &mesh_ref](const int in) {
        auto const vertex = mesh_ref.GetVertex(in);
        for (int d = 0; d < sdim; ++d) { nodeCoords_h(in * sdim + d) = vertex[d]; }
      });
      Kokkos::deep_copy(data.nodeCoords, nodeCoords_h);
    }

    // Copy cell-to-node connectivity
    data.cellToNode = Kokkos::View<int**, ExecutionSpace>("cellToNode", numCells, 4);  // Q1 has 4 nodes max
    {
      auto  cellToNode_h = Kokkos::create_mirror_view(data.cellToNode);
      auto& mesh_ref     = meshInfo.mesh;
      Kokkos::parallel_for("CellToNodes_Copy", Kokkos::RangePolicy<HostSpace>(0, numCells), [=, &mesh_ref](const int ic) {
        auto const& nodeList = mesh_ref.NodeList(ic);
        auto        c2n_ic   = Kokkos::subview(cellToNode_h, ic, Kokkos::ALL());
        for (int in = 0; in < nodeList.size() && in < 4; ++in) { c2n_ic(in) = nodeList[in]; }
      });
      Kokkos::deep_copy(data.cellToNode, cellToNode_h);
    }

    return data;
  }

 protected:
  /// \brief Host-side Q1 element assembly for MFEM fine elements
  ///
  /// Simplified version of ElementaryDataLagrangeFE_t for host execution
  /// Used by MFEM to assemble fine mesh elements
  ///
  void
  ElementaryDataLagrangeFE_Host(
      const double* __restrict__ coords_v,  // [8 values: x0,y0,x1,y1,x2,y2,x3,y3]
      double* __restrict__ rele,            // [4 values]
      double* __restrict__ kele) const      // [16 values]
  {
    constexpr int dim    = 2;
    constexpr int nNodes = fe2DQ1::numNode;

    // Initialize output
    for (int i = 0; i < nNodes; ++i) { rele[i] = 0.0; }
    for (int i = 0; i < nNodes * nNodes; ++i) { kele[i] = 0.0; }

    // Quadrature loop
    for (int iq = 0; iq < ruleLength; ++iq) {
      fe2DQ1 element;
      auto   NandGradN = element.GetValuesGradients(xi[iq], eta[iq], zeta[iq]);

      // Compute Jacobian
      double pointJac[dim * (dim + 1)] = {0};
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          double jacEntry = 0.0;
          for (int kn = 0; kn < nNodes; ++kn) { jacEntry += NandGradN[kn + jd * nNodes] * coords_v[id + kn * dim]; }
          pointJac[id + jd * dim] = jacEntry;
        }
      }

      double const xq = pointJac[0];
      double const yq = pointJac[1];

      // Evaluate material coefficients at quadrature point
      double alpha[dim];
      alpha[0] = ax_func(xq, yq, 0.0);
      alpha[1] = ay_func(xq, yq, 0.0);

      // Inverse Jacobian
      double detJ          = 1.0;
      double* __restrict J = &pointJac[dim];
      InverseInPlace<dim>(J, detJ);

      // Transform gradients
      double GradPhi[nNodes * dim];
      double* __restrict GradN = &NandGradN[nNodes];
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          double tmpGrad = 0.0;
          for (int kn = 0; kn < dim; ++kn) { tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
          GradPhi[in + jn * dim] = tmpGrad;
        }
      }

      // Assemble stiffness
      double w_v   = weight[iq];
      double coeff = w_v * detJ;

      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          double sum = 0.0;
          for (int kn = 0; kn < dim; ++kn) { sum += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim]; }
          kele[in + jn * nNodes] += sum * coeff;
        }
      }

      // Assemble RHS
      double fq = f_func(xq, yq, 0.0);
      for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * coeff; }
    }

    // Symmetrize the matrix once after all quadrature points
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
    }
  }

  /// \brief Host-side MFEM element assembly
  ///
  /// Implements multiscale FEM with static condensation
  /// Based on ElementaryDataMFEM_t from ScaledLaplacian.h
  ///
  void
  ElementaryDataMFEM_Host(
      const double*        coords_v,  // Coarse element coordinates [8 values]
      double*              rele,      // Coarse element RHS [4 values]
      double*              kele,      // Coarse element stiffness [16 values]
      std::vector<double>& phi)       // Basis functions (output)
  {
    constexpr int maxNumDofsPerEle = 4;
    int const     numNodes         = (ratio + 1) * (ratio + 1);

    std::vector<double> rhs(numNodes, 0);
    std::vector<double> rFineEle(maxNumDofsPerEle, 0);
    std::vector<double> kFineEle(maxNumDofsPerEle * maxNumDofsPerEle, 0);

    double const hx = (coords_v[2] - coords_v[0]) / double(ratio);
    double const hy = (coords_v[7] - coords_v[1]) / double(ratio);

    // Build fine mesh sparsity pattern
    std::vector<int> matRowPtr(numNodes + 1, 0);
    std::vector<int> matColIdx;
    matColIdx.reserve(9 * numNodes);

    for (int iy = 0; iy <= ratio; ++iy) {
      for (int ix = 0; ix <= ratio; ++ix) {
        int const nodeID = ix + iy * (ratio + 1);
        if (iy > 0) {
          if (ix > 0) { matColIdx.push_back(nodeID - 1 - (ratio + 1)); }
          matColIdx.push_back(nodeID - (ratio + 1));
          if (ix < ratio) { matColIdx.push_back(nodeID + 1 - (ratio + 1)); }
        }
        if (ix > 0) { matColIdx.push_back(nodeID - 1); }
        matColIdx.push_back(nodeID);
        if (ix < ratio) { matColIdx.push_back(nodeID + 1); }
        if (iy < ratio) {
          if (ix > 0) { matColIdx.push_back(nodeID - 1 + (ratio + 1)); }
          matColIdx.push_back(nodeID + (ratio + 1));
          if (ix < ratio) { matColIdx.push_back(nodeID + 1 + (ratio + 1)); }
        }
        matRowPtr[nodeID + 1] = matColIdx.size();
      }
    }

    std::vector<double> matValues(matColIdx.size(), 0);

    // Assemble fine elements (scalar path only, no SIMD)
    double coords[fe2DQ1::numNode * 2];
    for (int iy = 0; iy < ratio; ++iy) {
      for (int ix = 0; ix < ratio; ++ix) {
        std::array<int, fe2DQ1::numNode> nodeList{
            ix + iy * (ratio + 1),
            ix + 1 + iy * (ratio + 1),
            ix + 1 + (iy + 1) * (ratio + 1),
            ix + (iy + 1) * (ratio + 1)};

        coords[0] = coords_v[0] + ix * hx;
        coords[1] = coords_v[1] + iy * hy;
        coords[2] = coords[0] + hx;
        coords[3] = coords[1];
        coords[4] = coords[2];
        coords[5] = coords[3] + hy;
        coords[6] = coords[0];
        coords[7] = coords[5];

        std::fill(rFineEle.begin(), rFineEle.end(), 0);
        std::fill(kFineEle.begin(), kFineEle.end(), 0);

        this->ElementaryDataLagrangeFE_Host(&coords[0], &rFineEle[0], &kFineEle[0]);

        // Scatter to local system
        for (int in = 0; in < nodeList.size(); ++in) { rhs[nodeList[in]] += rFineEle[in]; }

        for (int in = 0; in < nodeList.size(); ++in) {
          auto const irow     = nodeList[in];
          auto const colBegin = &matColIdx[matRowPtr[irow]];
          auto const colEnd   = &matColIdx[matRowPtr[irow + 1]];
          for (int jn = 0; jn < nodeList.size(); ++jn) {
            auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
            matValues[matRowPtr[irow] + pos] += kFineEle[in + jn * nodeList.size()];
          }
        }
      }
    }

    // Create DOF mapping (interior vs boundary)
    std::vector<int> globalToFree(numNodes, -1);
    std::vector<int> freeToGlobal(numNodes - 4 * ratio);
    int              localCount = 0;
    for (int iy = 0; iy <= ratio; ++iy) {
      for (int ix = 0; ix <= ratio; ++ix) {
        if ((ix == 0) || (ix == ratio) || (iy == 0) || (iy == ratio)) { continue; }
        int const nodeID         = ix + iy * (ratio + 1);
        globalToFree[nodeID]     = localCount;
        freeToGlobal[localCount] = nodeID;
        localCount += 1;
      }
    }

    // Setup basis functions on boundaries
    int const numVectors = 4;
    phi.resize(numVectors * numNodes, 0);

    for (int iy = 0; iy <= ratio; ++iy) {
      double eta             = double(iy) / double(ratio);
      int    ix              = 0;
      int    in              = ix + iy * (ratio + 1);
      phi[in]                = 1.0 - eta;
      phi[in + 3 * numNodes] = eta;

      ix                     = ratio;
      in                     = ix + iy * (ratio + 1);
      phi[in + numNodes]     = 1.0 - eta;
      phi[in + 2 * numNodes] = eta;
    }

    for (int ix = 0; ix <= ratio; ++ix) {
      double xi          = double(ix) / double(ratio);
      int    iy          = 0;
      int    in          = ix + iy * (ratio + 1);
      phi[in]            = 1.0 - xi;
      phi[in + numNodes] = xi;

      iy                     = ratio;
      in                     = ix + iy * (ratio + 1);
      phi[in + 3 * numNodes] = 1.0 - xi;
      phi[in + 2 * numNodes] = xi;
    }

    // Build reduced system (without boundary DOFs)
    auto const       n = freeToGlobal.size();
    std::vector<int> newRowPtr(n + 1);
    newRowPtr[0] = 0;

    for (int i = 0; i < n; ++i) {
      auto gDof   = freeToGlobal[i];
      int  iCount = 0;
      for (auto k = matRowPtr[gDof]; k < matRowPtr[gDof + 1]; ++k) { iCount += (globalToFree[matColIdx[k]] != -1); }
      newRowPtr[i + 1] = newRowPtr[i] + iCount;
    }

    auto const          newNNZ = newRowPtr[n];
    std::vector<int>    newColIdx(newNNZ);
    std::vector<double> newValues(newNNZ);

    for (int iFree = 0; iFree < n; ++iFree) {
      auto const gdof = freeToGlobal[iFree];
      size_t     pos  = newRowPtr[iFree];
      for (auto k = matRowPtr[gdof]; k < matRowPtr[gdof + 1]; ++k) {
        auto const gCol = matColIdx[k];
        if (globalToFree[gCol] != -1) {
          newColIdx[pos] = globalToFree[gCol];
          newValues[pos] = matValues[k];
          pos += 1;
        }
      }
    }

    // Compute basis functions via static condensation
    std::vector<double>  btmp(numVectors * n, 0);
    SparseMatrix<double> K(numNodes, numNodes, matColIdx.size(), matRowPtr.data(), matColIdx.data(), matValues.data());
    std::vector<double>  Kphi(numVectors * numNodes);

    K.Apply(numVectors, &(phi[0]), &Kphi[0]);
    for (int ii = 0; ii < n; ++ii) {
      for (int ir = 0; ir < numVectors; ++ir) { btmp[ii + ir * n] = -Kphi[freeToGlobal[ii] + ir * numNodes]; }
    }

    // Solve for interior basis functions
    SymmetricSparse<double> Ktmp(n, newNNZ, newRowPtr.data(), newColIdx.data(), newValues.data());
    Ktmp.factor();

    std::vector<double> utmp(numVectors * n, 0);
    Ktmp.Solve(numVectors, &btmp[0], &utmp[0]);

    for (int ii = 0; ii < n; ++ii) {
      int grow = freeToGlobal[ii];
      for (int ir = 0; ir < numVectors; ++ir) { phi[grow + ir * numNodes] = utmp[ii + ir * n]; }
    }

    // Compute coarse element RHS
    for (int ir = 0; ir < numVectors; ++ir) {
      double sum = 0.0;
      for (int ii = 0; ii < numNodes; ++ii) { sum += phi[ii + ir * numNodes] * rhs[ii]; }
      rele[ir] = sum;
    }

    // Compute coarse element stiffness
    K.Apply(numVectors, &(phi[0]), &Kphi[0]);
    for (int ir = 0; ir < numVectors; ++ir) {
      for (int jr = 0; jr <= ir; ++jr) {
        double sum = 0.0;
        for (int ii = 0; ii < numNodes; ++ii) { sum += phi[ii + ir * numNodes] * Kphi[ii + jr * numNodes]; }
        kele[ir + jr * numVectors] = sum;
        kele[jr + ir * numVectors] = sum;
      }
    }
  }

  /// \brief Device kernel for Q1 element assembly
  ///
  /// Computes element stiffness matrix and RHS for a single Q1 element
  /// This function is designed to be called from within a CUDA kernel
  /// \param[in] coords Element nodal coordinates [8 values: x0,y0,x1,y1,x2,y2,x3,y3]
  /// \param[in] ax, ay, f Coefficient functors (types from class template)
  ///
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION static void
  ElementaryDataQ1(
      const Scalar* __restrict__ coords,
      const Scalar* __restrict__ quadWeight,
      const Scalar* __restrict__ quadXi,
      const Scalar* __restrict__ quadEta,
      const Scalar* __restrict__ quadZeta,
      int   ruleLen,
      FuncX ax,
      FuncY ay,
      FuncF f,
      Scalar* __restrict__ rele,  // Element RHS [4 values]
      Scalar* __restrict__ kele   // Element stiffness [16 values]
  )
  {
    constexpr int dim    = 2;
    constexpr int nNodes = fe2DQ1Cuda::numNode;

    // Initialize output arrays
    for (int i = 0; i < nNodes; ++i) { rele[i] = Scalar(0); }
    for (int i = 0; i < nNodes * nNodes; ++i) { kele[i] = Scalar(0); }

    Scalar NandGradN[nNodes * (dim + 1)];
    Scalar pointJac[dim * (dim + 1)];
    Scalar alpha[dim];
    Scalar GradPhi[nNodes * dim];

    // Quadrature loop
    for (int iq = 0; iq < ruleLen; ++iq) {
      fe2DQ1Cuda::GetValuesGradients(quadXi[iq], quadEta[iq], quadZeta[iq], NandGradN);

      // Compute Jacobian: pointJac = [x, y, dx/dxi, dy/dxi, dx/deta, dy/deta]
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          Scalar jacEntry = Scalar(0);
          for (int kn = 0; kn < nNodes; ++kn) { jacEntry += NandGradN[kn + jd * nNodes] * coords[id + kn * dim]; }
          pointJac[id + jd * dim] = jacEntry;
        }
      }

      auto const xq = pointJac[0];
      auto const yq = pointJac[1];

      // Get material coefficients (note: these are evaluated on host)
      alpha[0] = ax(xq, yq, 0);
      alpha[1] = ay(xq, yq, 0);

      // Compute inverse Jacobian
      Scalar detJ          = Scalar(1);
      Scalar* __restrict J = &pointJac[dim];
      InverseInPlaceCuda<dim>(J, detJ);

      // Transform gradients: GradPhi = J^T * GradN
      Scalar const* __restrict GradN = &NandGradN[nNodes];
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          Scalar tmpGrad = 0;
          for (int kn = 0; kn < dim; ++kn) { tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
          GradPhi[in + jn * dim] = tmpGrad;
        }
      }

      // Assemble element stiffness matrix (symmetric)
      Scalar w_v   = quadWeight[iq];
      Scalar coeff = w_v * detJ;
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          Scalar sum = Scalar(0);
          for (int kn = 0; kn < dim; ++kn) { sum += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim]; }
          kele[in + jn * nNodes] += sum * coeff;
        }
      }

      // Assemble RHS
      Scalar fq = f(xq, yq, 0);
      for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * coeff; }
    }  // for (int iq = 0; iq < ruleLen; ++iq)

    // Symmetrize the matrix once after all quadrature points
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
    }
  }
};

/// Implementation of GetLinearSystem
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<ExecutionSpace, FuncX, FuncY, FuncF>::GetLinearSystem(
    Kokkos::View<double*, ExecutionSpace> rhs,
    Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
    Kokkos::View<int*, ExecutionSpace>    matColIdx,
    Kokkos::View<double*, ExecutionSpace> matValues)
{
  // Get mesh info
  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) { throw std::runtime_error("Only 2D supported in CUDA version"); }

  printf("Number of colors: %d\n", c2e.numRows());

  // Copy mesh data to device
  auto meshData = CopyMeshToDevice();

  // Capture quadrature data for device
  auto quadWeight_local = quadWeight_d;
  auto quadXi_local     = quadXi_d;
  auto quadEta_local    = quadEta_d;
  auto quadZeta_local   = quadZeta_d;
  int  ruleLen          = ruleLength;

  // Capture coefficient functors for device (copy by value)
  auto ax_device = ax_func;
  auto ay_device = ay_func;
  auto f_device  = f_func;

  // Process each color
  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const numEle = c2e.row_map(ic + 1) - c2e.row_map(ic);
    if (numEle == 0) continue;
    printf("  Color %d: %d elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, ExecutionSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes      = meshData.cellTypes;
    auto nodeCoords     = meshData.nodeCoords;
    auto cellToNode     = meshData.cellToNode;
    auto eleList_device = eleList_d;

    // Parallel assembly for this color (no conflicts within same color)
    Kokkos::parallel_for(
        "Q1Assembly_Color", Kokkos::RangePolicy<ExecutionSpace>(0, numEle), KOKKOS_LAMBDA(const int ik) {
          auto const eleID = eleList_device(ik);

          // Only handle Q1 elements for now
          auto const cellType = static_cast<ElementType>(cellTypes(eleID));
          if (cellType != ElementType::Q1) {
            return;  // Skip non-Q1 elements
          }

          constexpr int nNodes = fe2DQ1Cuda::numNode;
          constexpr int dim    = 2;

          // Get locally element node indices and coordinates
          int    nodeList[nNodes];
          double coords[nNodes * dim];
          for (int i = 0; i < nNodes; ++i) {
            nodeList[i]         = cellToNode(eleID, i);
            coords[i * dim + 0] = nodeCoords(nodeList[i] * dim + 0);
            coords[i * dim + 1] = nodeCoords(nodeList[i] * dim + 1);
          }

          // Allocate element arrays
          double rele[nNodes];
          double kele[nNodes * nNodes];

          ElementaryDataQ1(
              &coords[0],
              &quadWeight_local(0),
              &quadXi_local(0),
              &quadEta_local(0),
              &quadZeta_local(0),
              ruleLen,
              ax_device,
              ay_device,
              f_device,
              &rele[0],
              &kele[0]);

          // Scatter to global arrays (no atomics needed due to coloring)
          for (int in = 0; in < nNodes; ++in) { rhs(nodeList[in]) += rele[in]; }

          for (int in = 0; in < nNodes; ++in) {
            auto const irow     = nodeList[in];
            auto const colBegin = &matColIdx(matRowPtr(irow));
            auto const colEnd   = &matColIdx(matRowPtr(irow + 1));
            for (int jn = 0; jn < nNodes; ++jn) {
              // Binary search for column position
              int const  numCols = colEnd - colBegin;
              auto const pos     = lower_bound_device(colBegin, numCols, nodeList[jn]);
              matValues(matRowPtr(irow) + pos) += kele[in + jn * nNodes];
            }
          }
        });

    Kokkos::fence();
  }
}

/// Implementation of GetLinearSystemMFEM
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<ExecutionSpace, FuncX, FuncY, FuncF>::GetLinearSystemMFEM(
    Kokkos::View<double*, ExecutionSpace> rhs,
    Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
    Kokkos::View<int*, ExecutionSpace>    matColIdx,
    Kokkos::View<double*, ExecutionSpace> matValues)
{
  // Get mesh info
  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) { throw std::runtime_error("Only 2D supported in CUDA version"); }

  printf("Number of colors: %d\n", c2e.numRows());

  // Copy mesh data to device
  auto meshData = CopyMeshToDevice();

  // Capture quadrature data for device
  auto quadWeight_local = quadWeight_d;
  auto quadXi_local     = quadXi_d;
  auto quadEta_local    = quadEta_d;
  auto quadZeta_local   = quadZeta_d;
  int  ruleLen          = ruleLength;

  // Capture coefficient functors for device (copy by value)
  auto ax_device = ax_func;
  auto ay_device = ay_func;
  auto f_device  = f_func;

  // Process each color
  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const numEle = c2e.row_map(ic + 1) - c2e.row_map(ic);
    if (numEle == 0) continue;
    printf("  Color %d: %d elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, ExecutionSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes      = meshData.cellTypes;
    auto nodeCoords     = meshData.nodeCoords;
    auto cellToNode     = meshData.cellToNode;
    auto eleList_device = eleList_d;

    // Calculate scratch memory requirements
    using functor_t = MFEMAssemblyFunctor<ExecutionSpace, FuncX, FuncY, FuncF>;
    using scratch_int_1d    = typename functor_t::scratch_int_1d;
    using scratch_double_1d = typename functor_t::scratch_double_1d;

    constexpr int numFineNodes = functor_t::numFineNodes;
    constexpr int maxNnzPerRow = 9;
    constexpr int maxNnz       = numFineNodes * maxNnzPerRow;

    // Scratch memory size calculation (level 1):
    // 2 double vectors: rhs, solution (for full fine system)
    // 4 int vectors for DOF mappings: globalToFree, freeToGlobal, globalToBoundary, boundaryToGlobal
    // 1 double matrix: phi (numVectors x numFineNodes)
    // 3 K_ii sparse matrix components: matRowPtr_ii, matColIdx_ii, matValues_ii
    // 3 K_b sparse matrix components: matRowPtr_b, matColIdx_b, matValues_b
    // 8 double vectors for static condensation: 4 RHS (btmp) + 4 solutions (utmp)
    // Note: We must allocate maxNnz_ii and maxNnz_b here since scratch size must be known at compile time.
    constexpr int numFreeDofs     = (functor_t::ratio - 1) * (functor_t::ratio - 1);
    constexpr int numBoundaryDofs = numFineNodes - numFreeDofs;  // 128 for ratio=32
    constexpr int numVectors      = 4;                           // Number of coarse element nodes
    constexpr int maxNnzPerRow_ii = 9;                           // Max neighbors for interior node (including diagonal)
    constexpr int maxNnz_ii       = numFreeDofs * maxNnzPerRow_ii;
    constexpr int maxNnzPerRow_b  = 9;                           // Max neighbors for boundary node
    constexpr int maxNnz_b        = numBoundaryDofs * maxNnzPerRow_b;
    size_t scratch_size_level1 =
        scratch_double_1d::shmem_size(numFineNodes) +                // rhs
        scratch_double_1d::shmem_size(numFineNodes) +                // solution
        scratch_int_1d::shmem_size(numFineNodes) +                   // globalToFree
        scratch_int_1d::shmem_size(numFreeDofs) +                    // freeToGlobal
        scratch_int_1d::shmem_size(numFineNodes) +                   // globalToBoundary
        scratch_int_1d::shmem_size(numBoundaryDofs) +                // boundaryToGlobal
        scratch_double_1d::shmem_size(numVectors * numFineNodes) +   // phi
        scratch_int_1d::shmem_size(numFreeDofs + 1) +                // matRowPtr_ii
        scratch_int_1d::shmem_size(maxNnz_ii) +                      // matColIdx_ii (conservative upper bound)
        scratch_double_1d::shmem_size(maxNnz_ii) +                   // matValues_ii (conservative upper bound)
        scratch_int_1d::shmem_size(numBoundaryDofs + 1) +            // matRowPtr_b
        scratch_int_1d::shmem_size(maxNnz_b) +                       // matColIdx_b (conservative upper bound)
        scratch_double_1d::shmem_size(maxNnz_b) +                    // matValues_b (conservative upper bound)
        scratch_double_1d::shmem_size(numFreeDofs * numVectors) +    // btmp (RHS vectors)
        scratch_double_1d::shmem_size(numFreeDofs * numVectors);     // utmp (solution vectors)

    // Create TeamPolicy with scratch memory
    Kokkos::TeamPolicy<ExecutionSpace> team_policy(numEle, Kokkos::AUTO);
    team_policy = team_policy.set_scratch_size(1, Kokkos::PerTeam(scratch_size_level1));

    // Parallel assembly for this color (no conflicts within same color)
    MFEMAssemblyFunctor<ExecutionSpace, FuncX, FuncY, FuncF> functor{
        eleList_device,
        cellTypes,
        nodeCoords,
        cellToNode,
        numEle,
        quadWeight_local,
        quadXi_local,
        quadEta_local,
        quadZeta_local,
        ruleLen,
        ax_device,
        ay_device,
        f_device};
    Kokkos::parallel_for("MFEMAssembly_Color", team_policy, functor);

    Kokkos::fence();
  }
}

/// Implementation of ProcessMFEMElements (host-side)
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<ExecutionSpace, FuncX, FuncY, FuncF>::ProcessMFEMElements(
    Kokkos::View<double*, HostSpace> rhs_h,
    Kokkos::View<size_t*, HostSpace> matRowPtr_h,
    Kokkos::View<int*, HostSpace>    matColIdx_h,
    Kokkos::View<double*, HostSpace> matValues_h)
{
  // Process MFEM_L elements sequentially on host
  // Uses ElementaryDataMFEM_Host for multiscale finite element assembly

  printf("Processing MFEM_L elements on host...\n");
  int mfemCount = 0;

  for (int eleID = 0; eleID < meshInfo.mesh.NumberCells(); ++eleID) {
    if (meshInfo.mesh.GetCellType(eleID) != ElementType::MFEM_L) { continue; }

    mfemCount++;
    auto const    nodeList = meshInfo.mesh.NodeList(eleID);
    constexpr int sdim     = 2;
    constexpr int nNodes   = fe2DQ1::numNode;

    // Get element coordinates
    double coords[nNodes * sdim];
    for (int i = 0; i < nNodes; ++i) {
      auto const vertex    = meshInfo.mesh.GetVertex(nodeList[i]);
      coords[i * sdim]     = vertex[0];
      coords[i * sdim + 1] = vertex[1];
    }

    // Allocate element arrays
    double rele[nNodes];
    double kele[nNodes * nNodes];
    for (int i = 0; i < nNodes; ++i) { rele[i] = 0.0; }
    for (int i = 0; i < nNodes * nNodes; ++i) { kele[i] = 0.0; }

    // Storage for basis functions
    std::vector<double> phi;

    // Call MFEM assembly
    printf("  Assembling MFEM_L element %d/%d...", mfemCount, eleID);
    fflush(stdout);
    this->ElementaryDataMFEM_Host(&coords[0], &rele[0], &kele[0], phi);
    printf(" done.\n");

    // Scatter to global RHS
    for (size_t in = 0; in < nodeList.size(); ++in) { rhs_h(nodeList[in]) += rele[in]; }

    // Scatter to global stiffness matrix
    for (size_t in = 0; in < nodeList.size(); ++in) {
      auto const irow     = nodeList[in];
      auto const colBegin = &matColIdx_h(matRowPtr_h(irow));
      auto const colEnd   = &matColIdx_h(matRowPtr_h(irow + 1));

      for (size_t jn = 0; jn < nodeList.size(); ++jn) {
        auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
        matValues_h(matRowPtr_h(irow) + pos) += kele[in + jn * nodeList.size()];
      }
    }
  }

  printf("Processed %d MFEM_L elements.\n", mfemCount);
}

}  // namespace IMSI
