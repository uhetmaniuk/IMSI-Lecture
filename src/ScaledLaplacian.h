#pragma once

#include <fstream>

#include "../main_config.h"
#include "Element.h"
#include "MathUtils.h"  // For InverseInPlace
#include "MathUtilsCuda.h"
#include "MeshUtils.h"
#include "PCG_Solver.h"
#include "QuadratureRule.h"
#include "SparseMatrix.hpp"
#include "fe2DQ1.h"

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
class ScaledLaplacian
{
 public:
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  ScaledLaplacian(
      const MeshConnectivity<>& meshData,
      RuleType                  quadRule,
      int                       quadOrder,
      FuncX                     ax_in,
      FuncY                     ay_in,
      FuncF                     f_in)
      : meshInfo(meshData),
        ruleType(quadRule),
        ruleOrder(quadOrder),
        ax_func(ax_in),
        ay_func(ay_in),
        f_func(f_in)
  {
    auto const sdim = meshInfo.mesh.GetSpatialDimension();
    if (sdim != 2) {
      throw std::runtime_error("ScaledLaplacian is only implemented for 2D problems");
    }

    std::vector<double> weight;
    std::vector<double> xi, eta, zeta;
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

  void
  OutputMFEMFine(const double* uCoarse, int numEleX, int numEleY) const;

  static constexpr int ratio = 32;  // For MFEM_L fine mesh refinement

 protected:
  const MeshConnectivity<> meshInfo;

  /// Coefficient functors (stored as member variables)
  FuncX ax_func;
  FuncY ay_func;
  FuncF f_func;

  RuleType ruleType   = RuleType::Gauss;
  int      ruleOrder  = 1;
  int      ruleLength = 0;

  // Device copies of quadrature data
  Kokkos::View<double*, ExecutionSpace> quadWeight_d;
  Kokkos::View<double*, ExecutionSpace> quadXi_d;
  Kokkos::View<double*, ExecutionSpace> quadEta_d;
  Kokkos::View<double*, ExecutionSpace> quadZeta_d;

  // MFEM basis functions stored on device [numElements, phiSize]
  // phiSize = numVectors * numFineNodes = 4 * (ratio+1)^2
  Kokkos::View<double**, ExecutionSpace> phiMFEM_d;

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
      Kokkos::parallel_for(
          "CellTypes_Copy", Kokkos::RangePolicy<HostSpace>(0, numCells), [=](const int ic) {
            cellTypes_h(ic) = static_cast<int>(cellTypes_ptr[ic]);
          });
      Kokkos::deep_copy(data.cellTypes, cellTypes_h);
    }

    // Copy node coordinates (interleaved)
    data.nodeCoords = Kokkos::View<double*, ExecutionSpace>("nodeCoords", numNodes * sdim);
    {
      auto  nodeCoords_h = Kokkos::create_mirror_view(data.nodeCoords);
      auto& mesh_ref     = meshInfo.mesh;
      Kokkos::parallel_for(
          "NodeCoords_Copy",
          Kokkos::RangePolicy<HostSpace>(0, numNodes),
          [=, &mesh_ref](const int in) {
            auto const vertex = mesh_ref.GetVertex(in);
            for (int d = 0; d < sdim; ++d) {
              nodeCoords_h(in * sdim + d) = vertex[d];
            }
          });
      Kokkos::deep_copy(data.nodeCoords, nodeCoords_h);
    }

    // Copy cell-to-node connectivity
    data.cellToNode =
        Kokkos::View<int**, ExecutionSpace>("cellToNode", numCells, 4);  // Q1 has 4 nodes max
    {
      auto  cellToNode_h = Kokkos::create_mirror_view(data.cellToNode);
      auto& mesh_ref     = meshInfo.mesh;
      Kokkos::parallel_for(
          "CellToNodes_Copy",
          Kokkos::RangePolicy<HostSpace>(0, numCells),
          [=, &mesh_ref](const int ic) {
            auto const& nodeList = mesh_ref.NodeList(ic);
            auto        c2n_ic   = Kokkos::subview(cellToNode_h, ic, Kokkos::ALL());
            for (int in = 0; in < nodeList.size() && in < 4; ++in) {
              c2n_ic(in) = nodeList[in];
            }
          });
      Kokkos::deep_copy(data.cellToNode, cellToNode_h);
    }

    return data;
  }

 public:
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
    /// TODO Could we generalize the routine by introducing ElementType as template parameter
    using ElementType = fe2DQ1;

    constexpr int dim    = 2;
    constexpr int nNodes = ElementType::numNode;

    // Initialize output arrays
    for (int i = 0; i < nNodes; ++i) {
      rele[i] = Scalar(0);
    }
    for (int i = 0; i < nNodes * nNodes; ++i) {
      kele[i] = Scalar(0);
    }

    Scalar NandGradN[nNodes * (dim + 1)];
    Scalar pointJac[dim * (dim + 1)];
    Scalar alpha[dim];
    Scalar GradPhi[nNodes * dim];

    // Quadrature loop
    for (int iq = 0; iq < ruleLen; ++iq) {
      ElementType::GetValuesGradients(quadXi[iq], quadEta[iq], quadZeta[iq], NandGradN);

      // Compute Jacobian: pointJac = [x, y, dx/dxi, dy/dxi, dx/deta, dy/deta]
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          Scalar jacEntry = Scalar(0);
          for (int kn = 0; kn < nNodes; ++kn) {
            jacEntry += NandGradN[kn + jd * nNodes] * coords[id + kn * dim];
          }
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
          for (int kn = 0; kn < dim; ++kn) {
            tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes];
          }
          GradPhi[in + jn * dim] = tmpGrad;
        }
      }

      // Assemble element stiffness matrix (symmetric)
      Scalar w_v   = quadWeight[iq];
      Scalar coeff = w_v * detJ;
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          Scalar sum = Scalar(0);
          for (int kn = 0; kn < dim; ++kn) {
            sum += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim];
          }
          kele[in + jn * nNodes] += sum * coeff;
        }
      }

      // Assemble RHS
      Scalar fq = f(xq, yq, 0);
      for (int in = 0; in < nNodes; ++in) {
        rele[in] += fq * NandGradN[in] * coeff;
      }
    }  // for (int iq = 0; iq < ruleLen; ++iq)

    // Symmetrize the matrix once after all quadrature points
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = jn + 1; in < nNodes; ++in) {
        kele[in + jn * nNodes] = kele[jn + in * nNodes];
      }
    }
  }
};

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
  size_t                                numEle;

  // Quadrature data
  Kokkos::View<double*, ExecutionSpace> quadWeight;
  Kokkos::View<double*, ExecutionSpace> quadXi;
  Kokkos::View<double*, ExecutionSpace> quadEta;
  Kokkos::View<double*, ExecutionSpace> quadZeta;
  int                                   ruleLen;

  // Global coarse-level RHS (output)
  Kokkos::View<double*, ExecutionSpace> globalRhs;

  // Global coarse-level stiffness matrix (CSR format, output)
  Kokkos::View<size_t*, ExecutionSpace> globalMatRowPtr;
  Kokkos::View<int*, ExecutionSpace>    globalMatColIdx;
  Kokkos::View<double*, ExecutionSpace> globalMatValues;

  // MFEM basis functions (output) [numElements, phiSize]
  Kokkos::View<double**, ExecutionSpace> phiMFEM_global;

  // Coefficient functors
  FuncX ax_func;
  FuncY ay_func;
  FuncF f_func;

  static int const ratio            = 32;
  static int constexpr numFineNodes = (ratio + 1) * (ratio + 1);
  static int constexpr numFineEle   = ratio * ratio;

  typedef typename Kokkos::TeamPolicy<ExecutionSpace>::member_type                   team_member;
  typedef typename ExecutionSpace::scratch_memory_space                              scratch_space;
  typedef Kokkos::View<int*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_1d;
  typedef Kokkos::View<double*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      scratch_double_1d;

  KOKKOS_FUNCTION
  void
  operator()(const team_member& teamMember) const
  {
    // Get element index
    auto const ieleCoarse = teamMember.league_rank();

    // Allocate scratch pad memory at team level (level 1)
    // 2 vectors of size numFineNodes
    scratch_double_1d rhs(teamMember.team_scratch(1), numFineNodes);

    // Extract globalToFree and freeToGlobal mappings for interior DOFs
    // Interior DOFs are those where (ix > 0) && (ix < ratio) && (iy > 0) && (iy < ratio)
    // Number of free DOFs = (ratio-1) * (ratio-1)
    auto const     numFreeDofs = (ratio - 1) * (ratio - 1);
    scratch_int_1d globalToFree(teamMember.team_scratch(1), numFineNodes);
    scratch_int_1d freeToGlobal(teamMember.team_scratch(1), numFreeDofs);

    // Extract globalToBoundary and boundaryToGlobal mappings for boundary DOFs
    // Boundary DOFs are those where (ix == 0) || (ix == ratio) || (iy == 0) || (iy == ratio)
    auto const     numBoundaryDofs = numFineNodes - numFreeDofs;  // 128 for ratio=32
    scratch_int_1d globalToBoundary(teamMember.team_scratch(1), numFineNodes);
    scratch_int_1d boundaryToGlobal(teamMember.team_scratch(1), numBoundaryDofs);

    // Setup basis functions on boundaries
    // phi is a matrix of size (numVectors x numFineNodes) stored in column-major order
    // Each column represents one basis function evaluated at all fine nodes
    constexpr int     numVectors = 4;  // Number of coarse element nodes
    scratch_double_1d phi(teamMember.team_scratch(1), numVectors * numFineNodes);

    // Initialize phi to zero
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, numVectors * numFineNodes), [&](int i) {
          phi(i) = 0;
        });

    // Initialize rhs to 0, globalToFree to -1 (indicating constrained/boundary DOF)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFineNodes), [&](int i) {
      rhs(i)              = 0;
      globalToFree(i)     = -1;
      globalToBoundary(i) = -1;
    });

    // 8 vectors for static condensation: 4 RHS + 4 solutions (one per basis function)
    scratch_double_1d btmp(teamMember.team_scratch(1), numFreeDofs * numVectors);  // RHS vectors
    // Solution vectors - no need to initialize utmp
    scratch_double_1d utmp(teamMember.team_scratch(1), numFreeDofs * numVectors);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs * numVectors), [&](int i) {
      btmp(i) = 0;
    });

    // Set a barrier as we will update phi, globalToFree, freeToGlobal, boundaryToGlobal,
    // globalToBoundary
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
            globalToFree(nodeID)    = freeCount;
            freeToGlobal(freeCount) = nodeID;
            freeCount += 1;
          } else {
            // Boundary node
            globalToBoundary(nodeID)        = boundaryCount;
            boundaryToGlobal(boundaryCount) = nodeID;
            boundaryCount += 1;
          }
        }
      }
      // Handle 4 corners for phi basis functions (avoids race conditions)
      // Corner (0, 0) - Basis 0
      phi(0 + 0 * numFineNodes) = 1.0;
      // Corner (ratio, 0) - Basis 1
      int in                     = ratio;
      phi(in + 1 * numFineNodes) = 1.0;
      // Corner (ratio, ratio) - Basis 2
      in                         = ratio + ratio * (ratio + 1);
      phi(in + 2 * numFineNodes) = 1.0;
      // Corner (0, ratio) - Basis 3
      in                         = ratio * (ratio + 1);
      phi(in + 3 * numFineNodes) = 1.0;
    });

    // Build boundary basis functions in parallel (skip corners)
    // Left and right edges (ix = 0 and ix = ratio)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 1, ratio), [&](int is) {
      auto const s = double(is) / double(ratio);  // abscissa along an edge
      // Left edge (ix = 0)
      int in                     = is * (ratio + 1);
      phi(in + 0 * numFineNodes) = 1.0 - s;  // Basis function 0
      phi(in + 3 * numFineNodes) = s;        // Basis function 3
      // Right edge (ix = ratio)
      in                         = ratio + is * (ratio + 1);
      phi(in + 1 * numFineNodes) = 1.0 - s;  // Basis function 1
      phi(in + 2 * numFineNodes) = s;        // Basis function 2
    });

    // Bottom and top edges (iy = 0 and iy = ratio, skip corners)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 1, ratio), [&](int is) {
      auto const s = double(is) / double(ratio);  // abscissa along an edge
      // Bottom edge (iy = 0)
      int in                     = is;
      phi(in + 0 * numFineNodes) = 1.0 - s;  // Basis function 0
      phi(in + 1 * numFineNodes) = s;        // Basis function 1
      // Top edge (iy = ratio)
      in                         = is + ratio * (ratio + 1);
      phi(in + 3 * numFineNodes) = 1.0 - s;  // Basis function 3
      phi(in + 2 * numFineNodes) = s;        // Basis function 2
    });

    // Build K_ii: CSR matrix for interior DOFs only
    // CSR matrix row pointer (we'll compute actual nnz before allocating colIdx and values)
    scratch_int_1d matRowPtr_ii(teamMember.team_scratch(1), numFreeDofs + 1);

    // Get the sparsity graph of K_ii (interior-interior stiffness matrix)
    // Map from free DOF index to sparsity pattern
    Kokkos::parallel_scan(
        Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int iFree, int& update, bool final) {
          int const iGlobal = freeToGlobal(iFree);
          int const ix      = iGlobal % (ratio + 1);
          int const iy      = iGlobal / (ratio + 1);

          // Count interior neighbors only
          int        count    = 1;  // Diagonal
          auto const hasWest  = (ix > 1);
          auto const hasEast  = (ix < ratio - 1);
          auto const hasSouth = (iy > 1);
          auto const hasNorth = (iy < ratio - 1);

          // South neighbors (iy - 1 row)
          if (hasSouth) {
            count += hasWest;  // SW
            count++;           // S (always present when hasSouth)
            count += hasEast;  // SE
          }
          // Center row neighbors
          count += hasWest + hasEast;  // W and E
          // North neighbors (iy + 1 row)
          if (hasNorth) {
            count += hasWest;  // NW
            count++;           // N (always present when hasNorth)
            count += hasEast;  // NE
          }

          if (final) {
            matRowPtr_ii(iFree) = update;
            if (iFree == numFreeDofs - 1)
              matRowPtr_ii(numFreeDofs) = update + count;
          }
          update += count;
        });
    teamMember.team_barrier();

    // Now allocate colIdx and values for K_ii with the exact number of non-zeros
    auto const        actualNnz_ii = matRowPtr_ii(numFreeDofs);
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
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
        // S neighbor
        {
          int jGlobal = iGlobal - (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
        // SE neighbor
        if (ix < ratio - 1) {
          int jGlobal = iGlobal + 1 - (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
      }
      // W neighbor
      if (ix > 1) {
        int jGlobal = iGlobal - 1;
        int jFree   = globalToFree(jGlobal);
        if (jFree != -1)
          matColIdx_ii(offset++) = jFree;
      }
      // Diagonal
      matColIdx_ii(offset++) = iFree;
      // E neighbor
      if (ix < ratio - 1) {
        int jGlobal = iGlobal + 1;
        int jFree   = globalToFree(jGlobal);
        if (jFree != -1)
          matColIdx_ii(offset++) = jFree;
      }
      // North neighbors (iy + 1 row)
      if (iy < ratio - 1) {
        // NW neighbor
        if (ix > 1) {
          int jGlobal = iGlobal - 1 + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
        // N neighbor
        {
          int jGlobal = iGlobal + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
        // NE neighbor
        if (ix < ratio - 1) {
          int jGlobal = iGlobal + 1 + (ratio + 1);
          int jFree   = globalToFree(jGlobal);
          if (jFree != -1)
            matColIdx_ii(offset++) = jFree;
        }
      }
    });

    // Initialize K_ii matrix values to zero
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, actualNnz_ii), [&](int i) {
      matValues_ii(i) = 0;
    });
    teamMember.team_barrier();

    // Build K_b: CSR matrix for boundary DOFs (numBoundaryDofs x numFineNodes)
    // K_b = [K_bi, K_bb] includes both boundary-interior and boundary-boundary couplings
    // We store K_b in global node numbering for columns (easier indexing)
    scratch_int_1d matRowPtr_b(teamMember.team_scratch(1), numBoundaryDofs + 1);

    // Get the sparsity graph of K_b (boundary rows, all columns)
    Kokkos::parallel_scan(
        Kokkos::TeamThreadRange(teamMember, numBoundaryDofs),
        [&](int iBoundary, int& update, bool final) {
          int const iGlobal = boundaryToGlobal(iBoundary);
          int const ix      = iGlobal % (ratio + 1);
          int const iy      = iGlobal / (ratio + 1);

          // Count all neighbors (interior and boundary) in global numbering
          int        count   = 1;  // Diagonal
          auto const hasWest = (ix > 0);
          auto const hasEast = (ix < ratio);
          count += hasWest + hasEast;  // W and E neighbors
          if (iy > 0) {
            count += 1 + hasWest + hasEast;
          }  // S, SW, SE neighbors
          if (iy < ratio) {
            count += 1 + hasWest + hasEast;
          }  // N, NW, NE neighbors

          if (final) {
            matRowPtr_b(iBoundary) = update;
            if (iBoundary == numBoundaryDofs - 1)
              matRowPtr_b(numBoundaryDofs) = update + count;
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
        if (ix > 0) {
          matColIdx_b(offset++) = iGlobal - 1 - (ratio + 1);
        }
        matColIdx_b(offset++) = iGlobal - (ratio + 1);
        if (ix < ratio) {
          matColIdx_b(offset++) = iGlobal + 1 - (ratio + 1);
        }
      }
      if (ix > 0) {
        matColIdx_b(offset++) = iGlobal - 1;
      }
      matColIdx_b(offset++) = iGlobal;  // Diagonal
      if (ix < ratio) {
        matColIdx_b(offset++) = iGlobal + 1;
      }
      if (iy < ratio) {
        if (ix > 0) {
          matColIdx_b(offset++) = iGlobal - 1 + (ratio + 1);
        }
        matColIdx_b(offset++) = iGlobal + (ratio + 1);
        if (ix < ratio) {
          matColIdx_b(offset++) = iGlobal + 1 + (ratio + 1);
        }
      }
    });

    // Initialize K_b matrix values to zero
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, actualNnz_b), [&](int i) {
      matValues_b(i) = 0;
    });
    teamMember.team_barrier();

    // Get coarse element coordinates from global mesh
    int const     eleID        = eleList_device(ieleCoarse);
    constexpr int nNodesCoarse = 4;  // Q1 coarse element
    constexpr int sdim         = 2;

    // Extract coarse element node coordinates
    double coords_coarse[nNodesCoarse * sdim];
    for (int i = 0; i < nNodesCoarse; ++i) {
      int const nodeID            = cellToNode(eleID, i);
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
      int           nodeList[nNodesFine];
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
      ScaledLaplacian<ExecutionSpace, FuncX, FuncY, FuncF>::template ElementaryDataQ1<double>(
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
              int const pos = lower_bound_device(&matColIdx_ii(rowBegin), rowEnd - rowBegin, jFree);
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
            int const pos = lower_bound_device(&matColIdx_b(rowBegin), rowEnd - rowBegin, jGlobal);
            Kokkos::atomic_add(&matValues_b(rowBegin + pos), k_val);
          }
        }
      }
    });

    teamMember.team_barrier();

    // ========================================================================
    // Solve the linear system using PCG with SSOR preconditioning
    // ========================================================================

    // Extract diagonal of K_ii for SSOR preconditioner
    scratch_double_1d diagValues_ii(teamMember.team_scratch(1), numFreeDofs);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int i) {
      diagValues_ii(i) = 0;
      for (int k = matRowPtr_ii(i); k < matRowPtr_ii(i + 1); ++k) {
        if (matColIdx_ii(k) == i) {
          diagValues_ii(i) = matValues_ii(k);
          break;
        }
      }
    });
    teamMember.team_barrier();

    // DEBUG: Print K_ii matrix info for element 0 only
    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
      if (ieleCoarse == 0) {
        printf("\n=== DEBUG: K_ii for element 0 ===\n");
        printf("numFreeDofs = %d\n", numFreeDofs);
        printf("actualNnz_ii = %d\n", matRowPtr_ii(numFreeDofs));

        // Print first few diagonal values
        printf("\nFirst 5 diagonal values:\n");
        for (int i = 0; i < 5 && i < numFreeDofs; ++i) {
          printf("  diag[%d] = %.6e\n", i, diagValues_ii(i));
        }

        // Print first few rows of K_ii
        printf("\nFirst 3 rows of K_ii (CSR):\n");
        for (int row = 0; row < 3 && row < numFreeDofs; ++row) {
          printf("  Row %d: ", row);
          for (int k = matRowPtr_ii(row); k < matRowPtr_ii(row + 1); ++k) {
            printf("(%d, %.4e) ", matColIdx_ii(k), matValues_ii(k));
          }
          printf("\n");
        }

        // Print first few btmp values (RHS for basis function 0)
        printf("\nFirst 5 btmp values (RHS for basis 0):\n");
        for (int i = 0; i < 5 && i < numFreeDofs; ++i) {
          printf("  btmp[%d] = %.6e\n", i, btmp(i));
        }

        // Check if matrix looks SPD (diagonal should be positive)
        int negDiag = 0;
        int zeroDiag = 0;
        for (int i = 0; i < numFreeDofs; ++i) {
          if (diagValues_ii(i) < 0) negDiag++;
          if (diagValues_ii(i) == 0) zeroDiag++;
        }
        printf("\nDiagonal check: %d negative, %d zero out of %d\n", negDiag, zeroDiag, numFreeDofs);
        printf("=== END DEBUG ===\n\n");
      }
    });
    teamMember.team_barrier();

    // Allocate workspace for PCG (4 * numFreeDofs per vector)
    constexpr int     numVectorsToSolve = 3;  // Solve only first 3, get 4th from partition of unity
    scratch_double_1d pcg_work(teamMember.team_scratch(1), 4 * numFreeDofs);

    // Initialize utmp with Q1 basis functions for good initial guess
    // For each interior node, compute parametric coordinates and evaluate Q1 basis
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int iFree) {
      int const  iGlobal   = freeToGlobal(iFree);
      int const  ix        = iGlobal % (ratio + 1);
      int const  iy        = iGlobal / (ratio + 1);
      auto const xi_param  = double(ix) / double(ratio);  // Parametric coords in [0,1]
      auto const eta_param = double(iy) / double(ratio);
      // Q1 basis functions: N0 = (1-xi)(1-eta), N1 = xi(1-eta), N2 = xi*eta
      utmp(iFree + 0 * numFreeDofs) = (1.0 - xi_param) * (1.0 - eta_param);  // Corner 0
      utmp(iFree + 1 * numFreeDofs) = xi_param * (1.0 - eta_param);          // Corner 1
      utmp(iFree + 2 * numFreeDofs) = xi_param * eta_param;                  // Corner 2
    });
    teamMember.team_barrier();

    // Solve for each of the first 3 basis functions using team-parallel PCG
    // All threads in the team participate in SpMV, dot products, and vector updates
    double const tol     = 1e-12;  // Convergence tolerance (tight for accuracy)
    int const    maxIter = 2000;   // Maximum PCG iterations (Jacobi needs more than SSOR)

    for (int ir = 0; ir < numVectorsToSolve; ++ir) {
      // Team-parallel PCG with Jacobi preconditioning
      // All threads participate - no Kokkos::single needed
      PCG_Solve_Team(
          teamMember,
          numFreeDofs,
          &matRowPtr_ii(0),
          &matColIdx_ii(0),
          &matValues_ii(0),
          &diagValues_ii(0),          // Diagonal for Jacobi preconditioner
          &btmp(ir * numFreeDofs),    // RHS
          &utmp(ir * numFreeDofs),    // Solution (initialized with Q1 basis)
          &pcg_work(0),               // Workspace (size 4*numFreeDofs)
          tol,
          maxIter);
      // Barrier after each solve to ensure completion before next
      teamMember.team_barrier();
    }

    // Reconstruct full basis functions by combining boundary (phi) and interior (utmp) values
    // Phi was initialized with boundary basis functions; now add interior solution
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numFreeDofs), [&](int iFree) {
      auto const iGlobal = freeToGlobal(iFree);
      double     sum     = 0;
      for (int ir = 0; ir < numVectorsToSolve; ++ir) {
        phi(iGlobal + ir * numFineNodes) = utmp(iFree + ir * numFreeDofs);
        sum += utmp(iFree + ir * numFreeDofs);
      }
      phi(iGlobal + 3 * numFineNodes) = 1.0 - sum;
    });
    teamMember.team_barrier();

    // Get coarse element ID for later scatter operations
    int const ieleCoarse_actual = eleList_device(ieleCoarse);

    // ========================================================================
    // Save phi basis functions to global memory for later reconstruction
    // ========================================================================
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, numVectors * numFineNodes), [&](int idx) {
          phiMFEM_global(ieleCoarse_actual, idx) = phi(idx);
        });
    teamMember.team_barrier();

    // ========================================================================
    // Compute coarse element stiffness matrix: kele = phi^T * Kfine * phi
    // ========================================================================

    // Allocate scratch memory for kele (4x4 symmetric matrix)
    constexpr int     nNodesCoarse_k = 4;  // Q1 coarse element
    scratch_double_1d kele(teamMember.team_scratch(1), nNodesCoarse_k * nNodesCoarse_k);

    // Initialize kele to zero
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, nNodesCoarse_k * nNodesCoarse_k), [&](int i) {
          kele(i) = 0;
        });
    teamMember.team_barrier();

    // Compute kele = phi^T * Kfine * phi
    //
    // kele = phi^T [[Kii phi_i + Kib phi_b]; [Kbi phi_i + Kbb phi_b]]
    // kele = phi^T [[0]; [Kbi phi_i + Kbb phi_b]]
    // kele = phi_b^T ( Kbi phi_i + Kbb phi_b )
    // kele = phi_b^T ( Kbb phi_b - Kbi Kii^{-1} Kib phi_b ) -> Schur complement

    // Contribution from K_b (boundary-all coupling)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, numBoundaryDofs), [&](int iBoundary) {
      int const iGlobal  = boundaryToGlobal(iBoundary);
      int const rowBegin = matRowPtr_b(iBoundary);
      int const rowEnd   = matRowPtr_b(iBoundary + 1);

      for (int k = rowBegin; k < rowEnd; ++k) {
        int const    jGlobal = matColIdx_b(k);  // K_b uses global node numbering for columns
        double const k_val   = matValues_b(k);

        // Accumulate contribution to kele for all basis function pairs
        for (int ir = 0; ir < numVectors; ++ir) {
          double const phi_i = phi(iGlobal + ir * numFineNodes);
          for (int jr = 0; jr < numVectors; ++jr) {
            double const phi_j = phi(jGlobal + jr * numFineNodes);
            Kokkos::atomic_add(&kele(ir + jr * nNodesCoarse_k), phi_i * k_val * phi_j);
          }
        }
      }
    });
    teamMember.team_barrier();

    // ========================================================================
    // Scatter-add kele into global coarse stiffness matrix
    // ========================================================================
    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
      for (int ir = 0; ir < numVectors; ++ir) {
        int const iCoarseNode = cellToNode(ieleCoarse_actual, ir);
        int const rowBegin    = globalMatRowPtr(iCoarseNode);
        int const rowEnd      = globalMatRowPtr(iCoarseNode + 1);

        for (int jr = 0; jr < numVectors; ++jr) {
          int const jCoarseNode = cellToNode(ieleCoarse_actual, jr);
          // Binary search for column position in global matrix
          int const pos =
              lower_bound_device(&globalMatColIdx(rowBegin), rowEnd - rowBegin, jCoarseNode);
          Kokkos::atomic_add(&globalMatValues(rowBegin + pos), kele(ir + jr * nNodesCoarse_k));
        }
      }
    });

    // ========================================================================
    // Compute coarse element RHS and scatter-add into global RHS
    // rele[ir] = phi[:,ir]^T * rhs_fine, then add to globalRhs
    // ========================================================================

    // Compute each component using team-level parallel reduction and scatter-add directly
    for (int ir = 0; ir < numVectors; ++ir) {
      double sum = 0;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(teamMember, numFineNodes),
          [&](int ii, double& local_sum) {
            local_sum += phi(ii + ir * numFineNodes) * rhs(ii);
          },
          sum);
      // Directly scatter-add into global coarse RHS
      Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
        int const coarseNode = cellToNode(ieleCoarse_actual, ir);
        Kokkos::atomic_add(&globalRhs(coarseNode), sum);
      });
    }
  }
};

/// Implementation of GetLinearSystem
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacian<ExecutionSpace, FuncX, FuncY, FuncF>::GetLinearSystem(
    Kokkos::View<double*, ExecutionSpace> rhs,
    Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
    Kokkos::View<int*, ExecutionSpace>    matColIdx,
    Kokkos::View<double*, ExecutionSpace> matValues)
{
  // Get mesh info
  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) {
    throw std::runtime_error("Only 2D supported in CUDA version");
  }

  printf("Number of colors: %ld\n", c2e.numRows());

  // Copy mesh data to device
  auto meshData_d = CopyMeshToDevice();

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
    if (numEle == 0)
      continue;
    printf("  Color %d: %ld elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, ExecutionSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(
          c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes      = meshData_d.cellTypes;
    auto nodeCoords     = meshData_d.nodeCoords;
    auto cellToNode     = meshData_d.cellToNode;
    auto eleList_device = eleList_d;

    // Parallel assembly for this color (no conflicts within same color)
    Kokkos::parallel_for(
        "Q1Assembly_Color",
        Kokkos::RangePolicy<ExecutionSpace>(0, numEle),
        KOKKOS_LAMBDA(const int ik) {
          auto const eleID = eleList_device(ik);

          // Only handle Q1 elements for now
          auto const cellType = static_cast<ElementType>(cellTypes(eleID));
          if (cellType != ElementType::Q1) {
            return;  // Skip non-Q1 elements
          }

          constexpr int nNodes = fe2DQ1::numNode;
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
          for (int in = 0; in < nNodes; ++in) {
            rhs(nodeList[in]) += rele[in];
          }

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
ScaledLaplacian<ExecutionSpace, FuncX, FuncY, FuncF>::GetLinearSystemMFEM(
    Kokkos::View<double*, ExecutionSpace> rhs,
    Kokkos::View<size_t*, ExecutionSpace> matRowPtr,
    Kokkos::View<int*, ExecutionSpace>    matColIdx,
    Kokkos::View<double*, ExecutionSpace> matValues)
{
  // Get mesh info
  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) {
    throw std::runtime_error("Only 2D supported in CUDA version");
  }

  printf("Number of colors: %ld\n", c2e.numRows());

  // Copy mesh data to device
  auto meshData_d = CopyMeshToDevice();

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

  // Allocate storage for MFEM basis functions on device
  constexpr int numVectors_const   = 4;  // Q1 coarse element has 4 nodes
  constexpr int numFineNodes_const = (ratio + 1) * (ratio + 1);
  constexpr int phiSize            = numVectors_const * numFineNodes_const;
  int const     numElements        = meshInfo.mesh.NumberCells();
  phiMFEM_d = Kokkos::View<double**, ExecutionSpace>("phiMFEM", numElements, phiSize);

  // Process each color
  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const numEle = c2e.row_map(ic + 1) - c2e.row_map(ic);
    if (numEle == 0)
      continue;
    printf("  Color %d: %ld elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, ExecutionSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(
          c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes      = meshData_d.cellTypes;
    auto nodeCoords     = meshData_d.nodeCoords;
    auto cellToNode     = meshData_d.cellToNode;
    auto eleList_device = eleList_d;

    // Calculate scratch memory requirements
    using functor_t         = MFEMAssemblyFunctor<ExecutionSpace, FuncX, FuncY, FuncF>;
    using scratch_int_1d    = typename functor_t::scratch_int_1d;
    using scratch_double_1d = typename functor_t::scratch_double_1d;

    constexpr int numFineNodes = functor_t::numFineNodes;

    // Scratch memory size calculation (level 1):
    // 2 double vectors: rhs, solution (for full fine system)
    // 4 int vectors for DOF mappings: globalToFree, freeToGlobal, globalToBoundary,
    // boundaryToGlobal 1 double matrix: phi (numVectors x numFineNodes) 3 K_ii sparse matrix
    // components: matRowPtr_ii, matColIdx_ii, matValues_ii 3 K_b sparse matrix components:
    // matRowPtr_b, matColIdx_b, matValues_b 8 double vectors for static condensation: 4 RHS (btmp)
    // + 4 solutions (utmp) 1 double vector: diagValues_ii (diagonal of K_ii for SSOR) 1 double
    // vector: pcg_work (workspace for PCG solver, size 4*numFreeDofs) 1 double matrix: kele (coarse
    // element stiffness matrix, 4x4 = 16 values) Note: We must allocate maxNnz_ii and maxNnz_b here
    // since scratch size must be known at compile time.
    constexpr int numFreeDofs     = (functor_t::ratio - 1) * (functor_t::ratio - 1);
    constexpr int numBoundaryDofs = numFineNodes - numFreeDofs;  // 128 for ratio=32
    constexpr int numVectors      = 4;                           // Number of coarse element nodes
    constexpr int maxNnzPerRow_ii = 9;  // Max neighbors for interior node (including diagonal)
    constexpr int maxNnz_ii       = numFreeDofs * maxNnzPerRow_ii;
    constexpr int maxNnzPerRow_b  = 9;  // Max neighbors for boundary node
    constexpr int maxNnz_b        = numBoundaryDofs * maxNnzPerRow_b;
    constexpr int nNodesCoarse    = 4;  // Q1 coarse element
    size_t        scratch_size_level1 =
        scratch_double_1d::shmem_size(numFineNodes) +               // rhs
        scratch_double_1d::shmem_size(numFineNodes) +               // solution
        scratch_int_1d::shmem_size(numFineNodes) +                  // globalToFree
        scratch_int_1d::shmem_size(numFreeDofs) +                   // freeToGlobal
        scratch_int_1d::shmem_size(numFineNodes) +                  // globalToBoundary
        scratch_int_1d::shmem_size(numBoundaryDofs) +               // boundaryToGlobal
        scratch_double_1d::shmem_size(numVectors * numFineNodes) +  // phi
        scratch_int_1d::shmem_size(numFreeDofs + 1) +               // matRowPtr_ii
        scratch_int_1d::shmem_size(maxNnz_ii) +     // matColIdx_ii (conservative upper bound)
        scratch_double_1d::shmem_size(maxNnz_ii) +  // matValues_ii (conservative upper bound)
        scratch_int_1d::shmem_size(numBoundaryDofs + 1) +  // matRowPtr_b
        scratch_int_1d::shmem_size(maxNnz_b) +             // matColIdx_b (conservative upper bound)
        scratch_double_1d::shmem_size(maxNnz_b) +          // matValues_b (conservative upper bound)
        scratch_double_1d::shmem_size(numFreeDofs * numVectors) +    // btmp (RHS vectors)
        scratch_double_1d::shmem_size(numFreeDofs * numVectors) +    // utmp (solution vectors)
        scratch_double_1d::shmem_size(numFreeDofs) +                 // diagValues_ii
        scratch_double_1d::shmem_size(4 * numFreeDofs) +             // pcg_work
        scratch_double_1d::shmem_size(nNodesCoarse * nNodesCoarse);  // kele

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
        rhs,
        matRowPtr,
        matColIdx,
        matValues,
        phiMFEM_d,
        ax_device,
        ay_device,
        f_device};
    Kokkos::parallel_for("MFEMAssembly_Color", team_policy, functor);

    Kokkos::fence();
  }
}

/// Implementation of OutputMFEMFine
template <typename ExecutionSpace, typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacian<ExecutionSpace, FuncX, FuncY, FuncF>::OutputMFEMFine(
    const double* uCoarse,
    int           numEleX,
    int           numEleY) const
{
  int const numFineEleX  = ratio * numEleX;
  int const numFineEleY  = ratio * numEleY;
  int const numFineNodes = (numFineEleX + 1) * (numFineEleY + 1);
  int const numCoarseEle = numEleX * numEleY;

  // Create Kokkos::Views for parallel computation on device
  Kokkos::View<double*, ExecutionSpace> uFine_d("uFine", numFineNodes);

  // Copy uCoarse from host to device (can't use unmanaged view with host pointer on GPU)
  int const numCoarseNodes = (numEleX + 1) * (numEleY + 1);
  Kokkos::View<double*, ExecutionSpace> uCoarse_d("uCoarse", numCoarseNodes);
  {
    Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        uCoarse_h(uCoarse, numCoarseNodes);
    Kokkos::deep_copy(uCoarse_d, uCoarse_h);
  }

  // Capture member variable for device lambda (avoid implicit 'this' capture)
  auto phiMFEM_local = phiMFEM_d;

  // Time the parallel reconstruction
  Kokkos::Timer timer;

  // Parallel reconstruction over coarse elements on device
  Kokkos::parallel_for(
      "MFEM_Reconstruction",
      Kokkos::RangePolicy<ExecutionSpace>(0, numCoarseEle),
      KOKKOS_LAMBDA(int eleID) {
        int const iy = eleID / numEleX;
        int const ix = eleID % numEleX;

        // Extract coarse DOFs for this element (uTrace)
        double uTrace[4];
        uTrace[0] = uCoarse_d(ix + iy * (numEleX + 1));
        uTrace[1] = uCoarse_d(ix + 1 + iy * (numEleX + 1));
        uTrace[2] = uCoarse_d(ix + 1 + (iy + 1) * (numEleX + 1));
        uTrace[3] = uCoarse_d(ix + (iy + 1) * (numEleX + 1));

        // Reconstruct fine nodes for this coarse element using stored phi
        // Precompute stride to avoid repeated multiplication
        constexpr int nodeStride = (ratio + 1) * (ratio + 1);
        for (int jy = 0; jy <= ratio; ++jy) {
          for (int jx = 0; jx <= ratio; ++jx) {
            int const localNode = jx + jy * (ratio + 1);
            int const nodeID    = ix * ratio + jx + (iy * ratio + jy) * (numFineEleX + 1);
            // Load all 4 phi values for this node (cache-friendly: consecutive access per k)
            double const phi0 = phiMFEM_local(eleID, localNode);
            double const phi1 = phiMFEM_local(eleID, localNode + nodeStride);
            double const phi2 = phiMFEM_local(eleID, localNode + 2 * nodeStride);
            double const phi3 = phiMFEM_local(eleID, localNode + 3 * nodeStride);
            uFine_d(nodeID) = phi0 * uTrace[0] + phi1 * uTrace[1] + phi2 * uTrace[2] + phi3 * uTrace[3];
          }
        }
      });
  Kokkos::fence();
  double computeTime = timer.seconds();

  // Copy result back to host for file output
  timer.reset();
  auto uFine_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), uFine_d);
  double copyTime = timer.seconds();

  // Write to file (use '\n' instead of std::endl to avoid flushing on each line)
  timer.reset();
  std::ofstream outFine("outputFine.txt");
  for (int i = 0; i < numFineNodes; ++i) {
    outFine << uFine_h(i) << '\n';
  }
  outFine.close();
  double ioTime = timer.seconds();

  std::cout << "    Compute time:    " << computeTime * 1000.0 << " ms\n";
  std::cout << "    D2H copy time:   " << copyTime * 1000.0 << " ms\n";
  std::cout << "    File I/O time:   " << ioTime * 1000.0 << " ms\n";
}

}  // namespace IMSI
