#pragma once

#include <functional>
#include <optional>

#include "../main_config.h"
#include "Element.h"
#include "MathUtilsCuda.h"
#include "MeshUtils.h"
#include "QuadratureRule.h"
#include "fe2DQ1Cuda.h"
#include "fe2DQ1.h"  // For host-side MFEM assembly
#include "MathUtils.h"  // For InverseInPlace
#include "SymmetricSparse.hpp"
#include "SparseMatrix.hpp"

namespace IMSI {

/// \brief Device-compatible binary search for sorted array
/// Returns index of found element or position where it would be inserted
template<typename T>
KOKKOS_INLINE_FUNCTION
int lower_bound_device(const T* array, int size, T value) {
  int left = 0;
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
///   - FuncX, FuncY, FuncF: Functor types for coefficients ax, ay, and f
///     Each must have operator()(double x, double y, double z) marked KOKKOS_INLINE_FUNCTION
///
template<typename FuncX, typename FuncY, typename FuncF>
class ScaledLaplacianCuda
{
 public:
  using CudaSpace = Kokkos::Cuda;
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  ScaledLaplacianCuda(
      const MeshConnectivity<>&                     meshData,
      RuleType                                      quadRule,
      int                                           quadOrder,
      FuncX                                         ax_in,
      FuncY                                         ay_in,
      FuncF                                         f_in)
      : meshInfo(meshData),
        ruleType(quadRule),
        ruleOrder(quadOrder),
        ax_func(ax_in),
        ay_func(ay_in),
        f_func(f_in)
  {
    auto const sdim = meshInfo.mesh.GetSpatialDimension();
    if (sdim != 2) {
      throw std::runtime_error("ScaledLaplacianCuda is only implemented for 2D problems");
    }

    getQuadrature(ruleType, sdim, ruleOrder, ruleLength, weight, xi, eta, zeta);

    // Copy quadrature data to device
    quadWeight_d = Kokkos::View<double*, CudaSpace>("quadWeight", ruleLength);
    quadXi_d = Kokkos::View<double*, CudaSpace>("quadXi", ruleLength);
    quadEta_d = Kokkos::View<double*, CudaSpace>("quadEta", ruleLength);
    quadZeta_d = Kokkos::View<double*, CudaSpace>("quadZeta", ruleLength);

    auto quadWeight_h = Kokkos::create_mirror_view(quadWeight_d);
    auto quadXi_h = Kokkos::create_mirror_view(quadXi_d);
    auto quadEta_h = Kokkos::create_mirror_view(quadEta_d);
    auto quadZeta_h = Kokkos::create_mirror_view(quadZeta_d);

    for (int i = 0; i < ruleLength; ++i) {
      quadWeight_h(i) = weight[i];
      quadXi_h(i) = xi[i];
      quadEta_h(i) = eta[i];
      quadZeta_h(i) = zeta[i];
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
      Kokkos::View<double*, CudaSpace> rhs,
      Kokkos::View<size_t*, CudaSpace> matRowPtr,
      Kokkos::View<int*, CudaSpace>    matColIdx,
      Kokkos::View<double*, CudaSpace> matValues);

  void
  GetLinearSystemMFEM(
      Kokkos::View<double*, CudaSpace> rhs,
      Kokkos::View<size_t*, CudaSpace> matRowPtr,
      Kokkos::View<int*, CudaSpace>    matColIdx,
      Kokkos::View<double*, CudaSpace> matValues);

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
  const MeshConnectivity<>                                     meshInfo;

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
  Kokkos::View<double*, CudaSpace> quadWeight_d;
  Kokkos::View<double*, CudaSpace> quadXi_d;
  Kokkos::View<double*, CudaSpace> quadEta_d;
  Kokkos::View<double*, CudaSpace> quadZeta_d;

 protected:
  /// \brief Host-side Q1 element assembly for MFEM fine elements
  ///
  /// Simplified version of ElementaryDataLagrangeFE_t for host execution
  /// Used by MFEM to assemble fine mesh elements
  ///
  void ElementaryDataLagrangeFE_Host(
      const double* __restrict__ coords_v,  // [8 values: x0,y0,x1,y1,x2,y2,x3,y3]
      double* __restrict__ rele,            // [4 values]
      double* __restrict__ kele) const      // [16 values]
  {
    constexpr int dim = 2;
    constexpr int nNodes = fe2DQ1::numNode;

    // Initialize output
    for (int i = 0; i < nNodes; ++i) { rele[i] = 0.0; }
    for (int i = 0; i < nNodes * nNodes; ++i) { kele[i] = 0.0; }

    // Quadrature loop
    for (int iq = 0; iq < ruleLength; ++iq) {
      fe2DQ1 element;
      auto NandGradN = element.GetValuesGradients(xi[iq], eta[iq], zeta[iq]);

      // Compute Jacobian
      double pointJac[dim * (dim + 1)] = {0};
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          double jacEntry = 0.0;
          for (int kn = 0; kn < nNodes; ++kn) {
            jacEntry += NandGradN[kn + jd * nNodes] * coords_v[id + kn * dim];
          }
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
      double detJ = 1.0;
      double* __restrict J = &pointJac[dim];
      InverseInPlace<dim>(J, detJ);

      // Transform gradients
      double GradPhi[nNodes * dim];
      double* __restrict GradN = &NandGradN[nNodes];
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          double tmpGrad = 0.0;
          for (int kn = 0; kn < dim; ++kn) {
            tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes];
          }
          GradPhi[in + jn * dim] = tmpGrad;
        }
      }

      // Assemble stiffness
      double w_v = weight[iq];
      double coeff = w_v * detJ;

      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          double sum = 0.0;
          for (int kn = 0; kn < dim; ++kn) {
            sum += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim];
          }
          kele[in + jn * nNodes] += sum * coeff;
        }
      }

      // Symmetrize
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = jn + 1; in < nNodes; ++in) {
          kele[in + jn * nNodes] = kele[jn + in * nNodes];
        }
      }

      // Assemble RHS
      double fq = f_func(xq, yq, 0.0);
      for (int in = 0; in < nNodes; ++in) {
        rele[in] += fq * NandGradN[in] * coeff;
      }
    }
  }

  /// \brief Host-side MFEM element assembly
  ///
  /// Implements multiscale FEM with static condensation
  /// Based on ElementaryDataMFEM_t from ScaledLaplacian.h
  ///
  void ElementaryDataMFEM_Host(
      const double* coords_v,     // Coarse element coordinates [8 values]
      double* rele,                // Coarse element RHS [4 values]
      double* kele,                // Coarse element stiffness [16 values]
      std::vector<double>& phi)    // Basis functions (output)
  {
    constexpr int maxNumDofsPerEle = 4;
    int const numNodes = (ratio + 1) * (ratio + 1);

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
        for (int in = 0; in < nodeList.size(); ++in) {
          rhs[nodeList[in]] += rFineEle[in];
        }

        for (int in = 0; in < nodeList.size(); ++in) {
          auto const irow = nodeList[in];
          auto const colBegin = &matColIdx[matRowPtr[irow]];
          auto const colEnd = &matColIdx[matRowPtr[irow + 1]];
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
    int localCount = 0;
    for (int iy = 0; iy <= ratio; ++iy) {
      for (int ix = 0; ix <= ratio; ++ix) {
        if ((ix == 0) || (ix == ratio) || (iy == 0) || (iy == ratio)) { continue; }
        int const nodeID = ix + iy * (ratio + 1);
        globalToFree[nodeID] = localCount;
        freeToGlobal[localCount] = nodeID;
        localCount += 1;
      }
    }

    // Setup basis functions on boundaries
    int const numVectors = 4;
    phi.resize(numVectors * numNodes, 0);

    for (int iy = 0; iy <= ratio; ++iy) {
      double eta = double(iy) / double(ratio);
      int ix = 0;
      int in = ix + iy * (ratio + 1);
      phi[in] = 1.0 - eta;
      phi[in + 3 * numNodes] = eta;

      ix = ratio;
      in = ix + iy * (ratio + 1);
      phi[in + numNodes] = 1.0 - eta;
      phi[in + 2 * numNodes] = eta;
    }

    for (int ix = 0; ix <= ratio; ++ix) {
      double xi = double(ix) / double(ratio);
      int iy = 0;
      int in = ix + iy * (ratio + 1);
      phi[in] = 1.0 - xi;
      phi[in + numNodes] = xi;

      iy = ratio;
      in = ix + iy * (ratio + 1);
      phi[in + 3 * numNodes] = 1.0 - xi;
      phi[in + 2 * numNodes] = xi;
    }

    // Build reduced system (without boundary DOFs)
    auto const n = freeToGlobal.size();
    std::vector<int> newRowPtr(n + 1);
    newRowPtr[0] = 0;

    for (int i = 0; i < n; ++i) {
      auto gDof = freeToGlobal[i];
      int iCount = 0;
      for (auto k = matRowPtr[gDof]; k < matRowPtr[gDof + 1]; ++k) {
        iCount += (globalToFree[matColIdx[k]] != -1);
      }
      newRowPtr[i + 1] = newRowPtr[i] + iCount;
    }

    auto const newNNZ = newRowPtr[n];
    std::vector<int> newColIdx(newNNZ);
    std::vector<double> newValues(newNNZ);

    for (int iFree = 0; iFree < n; ++iFree) {
      auto const gdof = freeToGlobal[iFree];
      size_t pos = newRowPtr[iFree];
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
    std::vector<double> btmp(numVectors * n, 0);
    SparseMatrix<double> K(
        numNodes, numNodes, matColIdx.size(), matRowPtr.data(), matColIdx.data(), matValues.data());
    std::vector<double> Kphi(numVectors * numNodes);

    K.Apply(numVectors, &(phi[0]), &Kphi[0]);
    for (int ii = 0; ii < n; ++ii) {
      for (int ir = 0; ir < numVectors; ++ir) {
        btmp[ii + ir * n] = -Kphi[freeToGlobal[ii] + ir * numNodes];
      }
    }

    // Solve for interior basis functions
    SymmetricSparse<double> Ktmp(n, newNNZ, newRowPtr.data(), newColIdx.data(), newValues.data());
    Ktmp.factor();

    std::vector<double> utmp(numVectors * n, 0);
    Ktmp.Solve(numVectors, &btmp[0], &utmp[0]);

    for (int ii = 0; ii < n; ++ii) {
      int grow = freeToGlobal[ii];
      for (int ir = 0; ir < numVectors; ++ir) {
        phi[grow + ir * numNodes] = utmp[ii + ir * n];
      }
    }

    // Compute coarse element RHS
    for (int ir = 0; ir < numVectors; ++ir) {
      double sum = 0.0;
      for (int ii = 0; ii < numNodes; ++ii) {
        sum += phi[ii + ir * numNodes] * rhs[ii];
      }
      rele[ir] = sum;
    }

    // Compute coarse element stiffness
    K.Apply(numVectors, &(phi[0]), &Kphi[0]);
    for (int ir = 0; ir < numVectors; ++ir) {
      for (int jr = 0; jr <= ir; ++jr) {
        double sum = 0.0;
        for (int ii = 0; ii < numNodes; ++ii) {
          sum += phi[ii + ir * numNodes] * Kphi[ii + jr * numNodes];
        }
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
  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION
  static void ElementaryDataQ1(
      const Scalar* __restrict__ coords,
      const Scalar* __restrict__  quadWeight,
      const Scalar* __restrict__  quadXi,
      const Scalar* __restrict__  quadEta,
      const Scalar* __restrict__  quadZeta,
      int ruleLen,
      FuncX ax, FuncY ay, FuncF f,
      Scalar* __restrict__ rele,           // Element RHS [4 values]
      Scalar* __restrict__ kele           // Element stiffness [16 values]
)
  {
    constexpr int dim = 2;
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
      Scalar detJ = Scalar(1);
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
      Scalar w_v = quadWeight[iq];
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

      // Symmetrize
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = jn + 1; in < nNodes; ++in) {
          kele[in + jn * nNodes] = kele[jn + in * nNodes];
        }
      }

      // Assemble RHS
      Scalar fq = f(xq, yq, 0);
      for (int in = 0; in < nNodes; ++in) {
        rele[in] += fq * NandGradN[in] * coeff;
      }
    } // for (int iq = 0; iq < ruleLen; ++iq)
  }
};

/// Implementation of GetLinearSystem
template<typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<FuncX, FuncY, FuncF>::GetLinearSystem(
    Kokkos::View<double*, CudaSpace> rhs,
    Kokkos::View<size_t*, CudaSpace> matRowPtr,
    Kokkos::View<int*, CudaSpace>    matColIdx,
    Kokkos::View<double*, CudaSpace> matValues)
{
  // Get mesh info
  auto const& c2e = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) {
    throw std::runtime_error("Only 2D supported in CUDA version");
  }

  printf("Number of colors: %d\n", c2e.numRows());

  // ====================================================================
  // Extract mesh data to device-accessible structures
  // ====================================================================

  int const numCells = meshInfo.mesh.NumberCells();
  int const numNodes = meshInfo.mesh.NumberVertices();

  // Create device views for mesh data
  Kokkos::View<int*, CudaSpace> cellTypes_d("cellTypes", numCells);
  {
    auto cellTypes_h = Kokkos::create_mirror_view(cellTypes_d);
    auto *cellTypes_ptr = meshInfo.mesh.GetCellType().data();
    Kokkos::parallel_for(
      "CellTypes_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numCells),
      [=](const int ic) {
        cellTypes_h(ic) = static_cast<int>(cellTypes_ptr[ic]);
      });
    Kokkos::deep_copy(cellTypes_d, cellTypes_h);
  }

  Kokkos::View<double*, CudaSpace> nodeCoords_d("nodeCoords", numNodes * sdim);  // x,y interleaved
  {
    auto nodeCoords_h = Kokkos::create_mirror_view(nodeCoords_d);
    auto& mesh_ref = meshInfo.mesh;
    Kokkos::parallel_for(
      "NodeCoords_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numNodes),
      [=, &mesh_ref](const int in) {
        auto const vertex = mesh_ref.GetVertex(in);
        nodeCoords_h(in * 2 + 0) = vertex[0];
        nodeCoords_h(in * 2 + 1) = vertex[1];
      });
    Kokkos::deep_copy(nodeCoords_d, nodeCoords_h);
  }

  Kokkos::View<int**, CudaSpace> cellToNode_d("cellToNode", numCells, 4); // Q1 has 4 nodes max
  {
    auto cellToNode_h = Kokkos::create_mirror_view(cellToNode_d);
    auto& mesh_ref = meshInfo.mesh;
    Kokkos::parallel_for(
      "CellToNodes_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numCells),
      [=, &mesh_ref](const int ic) {
        auto const& nodeList = mesh_ref.NodeList(ic);
        auto c2n_ic = Kokkos::subview(cellToNode_h, ic, Kokkos::ALL());
        for (int in = 0; in < nodeList.size() && in < 4; ++in) {
          c2n_ic(in) = nodeList[in];
        }
      });
    Kokkos::deep_copy(cellToNode_d, cellToNode_h);
  }

  // Capture quadrature data for device
  auto quadWeight_local = quadWeight_d;
  auto quadXi_local = quadXi_d;
  auto quadEta_local = quadEta_d;
  auto quadZeta_local = quadZeta_d;
  int ruleLen = ruleLength;

  // Capture coefficient functors for device (copy by value)
  auto ax_device = ax_func;
  auto ay_device = ay_func;
  auto f_device = f_func;

  // Process each color
  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const numEle = c2e.row_map(ic + 1) - c2e.row_map(ic);
    if (numEle == 0) continue;
    printf("  Color %d: %d elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, CudaSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes = cellTypes_d;
    auto nodeCoords = nodeCoords_d;
    auto cellToNode = cellToNode_d;
    auto eleList_device = eleList_d;

    // Parallel assembly for this color (no conflicts within same color)
    Kokkos::parallel_for(
        "Q1Assembly_Color",
        Kokkos::RangePolicy<CudaSpace>(0, numEle),
        KOKKOS_LAMBDA(const int ik) {
          auto const eleID = eleList_device(ik);

          // Only handle Q1 elements for now
          auto const cellType = static_cast<ElementType>(cellTypes(eleID));
          if (cellType != ElementType::Q1) {
            return;  // Skip non-Q1 elements
          }

          constexpr int nNodes = fe2DQ1Cuda::numNode;
          constexpr int dim = 2;

          // Get locally element node indices and coordinates
          int nodeList[nNodes];
          double coords[nNodes * dim];
          for (int i = 0; i < nNodes; ++i) {
            nodeList[i] = cellToNode(eleID, i);
            coords[i * dim + 0] = nodeCoords(nodeList[i] * dim + 0);
            coords[i * dim + 1] = nodeCoords(nodeList[i] * dim + 1);
          }

          // Allocate element arrays
          double rele[nNodes];
          double kele[nNodes * nNodes];

          ElementaryDataQ1(&coords[0], &quadWeight_local(0), &quadXi_local(0), &quadEta_local(0), &quadZeta_local(0), ruleLen,
            ax_device, ay_device, f_device, &rele[0], &kele[0]);

      /*
          // Compute element matrices
          for (int i = 0; i < nNodes; ++i) { rele[i] = 0.0; }
          for (int i = 0; i < nNodes * nNodes; ++i) { kele[i] = 0.0; }

          // === Inline Q1 assembly ===
          double NandGradN[nNodes * (dim + 1)];
          double pointJac[dim * (dim + 1)];
          double GradPhi[nNodes * dim];

          for (int iq = 0; iq < ruleLen; ++iq) {
            fe2DQ1Cuda::GetValuesGradients(
                quadXi_local(iq), quadEta_local(iq), quadZeta_local(iq), NandGradN);

            // Compute Jacobian
            for (int jd = 0; jd <= dim; ++jd) {
              for (int id = 0; id < dim; ++id) {
                double jacEntry = 0.0;
                for (int kn = 0; kn < nNodes; ++kn) {
                  jacEntry += NandGradN[kn + jd * nNodes] * coords[id + kn * dim];
                }
                pointJac[id + jd * dim] = jacEntry;
              }
            }

            // Extract physical coordinates of quadrature point
            double const xq = pointJac[0];
            double const yq = pointJac[1];

            // Inverse Jacobian
            double detJ = 1.0;
            double* __restrict J = &pointJac[dim];
            InverseInPlaceCuda<dim>(J, detJ);

            // Transform gradients
            double* __restrict GradN = &NandGradN[nNodes];
            for (int jn = 0; jn < nNodes; ++jn) {
              for (int in = 0; in < dim; ++in) {
                double tmpGrad = 0.0;
                for (int kn = 0; kn < dim; ++kn) {
                  tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes];
                }
                GradPhi[in + jn * dim] = tmpGrad;
              }
            }

            // Evaluate material coefficients at quadrature point (device-side)
            // TODO: Make this more general - currently hardcoded for specific functions
            double alpha_x_q = 1.0;  // Constant for now
            double alpha_y_q = 1.0;  // Constant for now
            // For varying coefficients, evaluate at (xq, yq):
            // double alpha_x_q = ... function of xq, yq ...
            // double alpha_y_q = ... function of xq, yq ...

            // Assemble stiffness
            double w_v = quadWeight_local(iq);
            double coeff = w_v * detJ;

            for (int jn = 0; jn < nNodes; ++jn) {
              auto const a_dphi_xj = GradPhi[0 + jn * dim] * alpha_x_q;
              auto const a_dphi_yj = GradPhi[1 + jn * dim] * alpha_y_q;
              for (int in = jn; in <= nNodes; ++in) {
                double sum = 0.0;
                sum += GradPhi[0 + in * dim] * a_dphi_xj;
                sum += GradPhi[1 + in * dim] * a_dphi_yj;
                kele[in + jn * nNodes] += sum * coeff;
              }
            }

            // Symmetrize
            for (int jn = 0; jn < nNodes; ++jn) {
              for (int in = 0; in < jn; ++in) {
                kele[in + jn * nNodes] = kele[jn + in * nNodes];
              }
            }

            // Assemble RHS
            if (has_f) {
              // Evaluate forcing function at quadrature point (device-side)
              // TODO: Make this more general - currently hardcoded for specific function
              // Current function: f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
              constexpr double pi = 3.14159265358979323846;
              double fq = 2.0 * pi * pi * sin(pi * xq) * sin(pi * yq);
              for (int in = 0; in < nNodes; ++in) {
                rele[in] += fq * NandGradN[in] * coeff;
              }
            }
          }
      */
          // === End inline Q1 assembly ===

          // Scatter to global arrays (no atomics needed due to coloring)
          for (int in = 0; in < nNodes; ++in) {
            rhs(nodeList[in]) += rele[in];
          }

          for (int in = 0; in < nNodes; ++in) {
            auto const irow = nodeList[in];
            auto const colBegin = &matColIdx(matRowPtr(irow));
            auto const colEnd = &matColIdx(matRowPtr(irow + 1));
            for (int jn = 0; jn < nNodes; ++jn) {
              // Binary search for column position
              int const numCols = colEnd - colBegin;
              auto const pos = lower_bound_device(colBegin, numCols, nodeList[jn]);
              matValues(matRowPtr(irow) + pos) += kele[in + jn * nNodes];
            }
          }
        });

    Kokkos::fence();
  }

  /*
  // Handle MFEM_L elements if present (process on host)
  bool hasMFEM = false;
  for (int ic = 0; ic < meshInfo.mesh.NumberCells(); ++ic) {
    if (meshInfo.mesh.GetCellType(ic) == ElementType::MFEM_L) {
      hasMFEM = true;
      break;
    }
  }

  if (hasMFEM) {
    printf("Processing MFEM_L elements on host...\n");

    // Transfer data to host
    auto rhs_h = Kokkos::create_mirror_view(rhs);
    auto matRowPtr_h = Kokkos::create_mirror_view(matRowPtr);
    auto matColIdx_h = Kokkos::create_mirror_view(matColIdx);
    auto matValues_h = Kokkos::create_mirror_view(matValues);

    Kokkos::deep_copy(rhs_h, rhs);
    Kokkos::deep_copy(matRowPtr_h, matRowPtr);
    Kokkos::deep_copy(matColIdx_h, matColIdx);
    Kokkos::deep_copy(matValues_h, matValues);

    // Process MFEM elements on host
    ProcessMFEMElements(rhs_h, matRowPtr_h, matColIdx_h, matValues_h);

    // Copy back to device
    Kokkos::deep_copy(rhs, rhs_h);
    Kokkos::deep_copy(matValues, matValues_h);

    printf("MFEM_L processing complete.\n");
  }
  */
}

/// Implementation of GetLinearSystem
template<typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<FuncX, FuncY, FuncF>::GetLinearSystemMFEM(
    Kokkos::View<double*, CudaSpace> rhs,
    Kokkos::View<size_t*, CudaSpace> matRowPtr,
    Kokkos::View<int*, CudaSpace>    matColIdx,
    Kokkos::View<double*, CudaSpace> matValues)
{
  // Get mesh info
  auto const& c2e = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  if (sdim != 2) {
    throw std::runtime_error("Only 2D supported in CUDA version");
  }

  printf("Number of colors: %d\n", c2e.numRows());

  // ====================================================================
  // Extract mesh data to device-accessible structures
  // ====================================================================

  int const numCells = meshInfo.mesh.NumberCells();
  int const numNodes = meshInfo.mesh.NumberVertices();

  // Create device views for mesh data
  Kokkos::View<int*, CudaSpace> cellTypes_d("cellTypes", numCells);
  {
    auto cellTypes_h = Kokkos::create_mirror_view(cellTypes_d);
    auto *cellTypes_ptr = meshInfo.mesh.GetCellType().data();
    Kokkos::parallel_for(
      "CellTypes_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numCells),
      [=](const int ic) {
        cellTypes_h(ic) = static_cast<int>(cellTypes_ptr[ic]);
      });
    Kokkos::deep_copy(cellTypes_d, cellTypes_h);
  }

  Kokkos::View<double*, CudaSpace> nodeCoords_d("nodeCoords", numNodes * sdim);  // x,y interleaved
  {
    auto nodeCoords_h = Kokkos::create_mirror_view(nodeCoords_d);
    auto& mesh_ref = meshInfo.mesh;
    Kokkos::parallel_for(
      "NodeCoords_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numNodes),
      [=, &mesh_ref](const int in) {
        auto const vertex = mesh_ref.GetVertex(in);
        nodeCoords_h(in * 2 + 0) = vertex[0];
        nodeCoords_h(in * 2 + 1) = vertex[1];
      });
    Kokkos::deep_copy(nodeCoords_d, nodeCoords_h);
  }

  Kokkos::View<int**, CudaSpace> cellToNode_d("cellToNode", numCells, 4); // Q1 has 4 nodes max
  {
    auto cellToNode_h = Kokkos::create_mirror_view(cellToNode_d);
    auto& mesh_ref = meshInfo.mesh;
    Kokkos::parallel_for(
      "CellToNodes_Copy",
      Kokkos::RangePolicy<HostSpace>(0, numCells),
      [=, &mesh_ref](const int ic) {
        auto const& nodeList = mesh_ref.NodeList(ic);
        auto c2n_ic = Kokkos::subview(cellToNode_h, ic, Kokkos::ALL());
        for (int in = 0; in < nodeList.size() && in < 4; ++in) {
          c2n_ic(in) = nodeList[in];
        }
      });
    Kokkos::deep_copy(cellToNode_d, cellToNode_h);
  }

  // Capture quadrature data for device
  auto quadWeight_local = quadWeight_d;
  auto quadXi_local = quadXi_d;
  auto quadEta_local = quadEta_d;
  auto quadZeta_local = quadZeta_d;
  int ruleLen = ruleLength;

  // Capture coefficient functors for device (copy by value)
  auto ax_device = ax_func;
  auto ay_device = ay_func;
  auto f_device = f_func;

  // Process each color
  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const numEle = c2e.row_map(ic + 1) - c2e.row_map(ic);
    if (numEle == 0) continue;
    printf("  Color %d: %d elements\n", ic, numEle);

    // Copy element list to device (eleList from rowConst is host-side!)
    Kokkos::View<int*, CudaSpace> eleList_d("eleList", numEle);
    {
      auto eleList_h = Kokkos::subview(c2e.entries, Kokkos::make_pair(c2e.row_map(ic), c2e.row_map(ic) + numEle));
      Kokkos::deep_copy(eleList_d, eleList_h);
    }

    // Capture device views for lambda
    auto cellTypes = cellTypes_d;
    auto nodeCoords = nodeCoords_d;
    auto cellToNode = cellToNode_d;
    auto eleList_device = eleList_d;

    int const warpSize = 32;
    int const teamSize = 32 * warpSize;
    int len = (numEle + teamSize - 1) / teamSize;
    Kokkos::TeamPolicy<CudaSpace> team_policy(len, teamSize, warpSize);
    // Parallel assembly for this color (no conflicts within same color)
    Kokkos::parallel_for(
        "MFEMAssembly_Color",
        team_policy,
        KOKKOS_LAMBDA(const Kokkos::TeamPolicy<CudaSpace>::member_type& teamMember) {
          const int ik = (team_member.league_rank() * team_member.team_size() + teamMember.team_rank()) / warpSize;
          if (ik < numEle) {
            auto const eleID = eleList_device(ik);
          }
          //
///          Kokkos::parallel_for("QQQQ", )
      /*
          // Scatter to global arrays (no atomics needed due to coloring)
          for (int in = 0; in < nNodes; ++in) {
            rhs(nodeList[in]) += rele[in];
          }

          for (int in = 0; in < nNodes; ++in) {
            auto const irow = nodeList[in];
            auto const colBegin = &matColIdx(matRowPtr(irow));
            auto const colEnd = &matColIdx(matRowPtr(irow + 1));
            for (int jn = 0; jn < nNodes; ++jn) {
              // Binary search for column position
              int const numCols = colEnd - colBegin;
              auto const pos = lower_bound_device(colBegin, numCols, nodeList[jn]);
              matValues(matRowPtr(irow) + pos) += kele[in + jn * nNodes];
            }
          }
          */

        });

    Kokkos::fence();
  }

}

/// Implementation of ProcessMFEMElements (host-side)
template<typename FuncX, typename FuncY, typename FuncF>
void
ScaledLaplacianCuda<FuncX, FuncY, FuncF>::ProcessMFEMElements(
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
    if (meshInfo.mesh.GetCellType(eleID) != ElementType::MFEM_L) {
      continue;
    }

    mfemCount++;
    auto const nodeList = meshInfo.mesh.NodeList(eleID);
    constexpr int sdim = 2;
    constexpr int nNodes = fe2DQ1::numNode;

    // Get element coordinates
    double coords[nNodes * sdim];
    for (int i = 0; i < nNodes; ++i) {
      auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
      coords[i * sdim] = vertex[0];
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
    for (size_t in = 0; in < nodeList.size(); ++in) {
      rhs_h(nodeList[in]) += rele[in];
    }

    // Scatter to global stiffness matrix
    for (size_t in = 0; in < nodeList.size(); ++in) {
      auto const irow = nodeList[in];
      auto const colBegin = &matColIdx_h(matRowPtr_h(irow));
      auto const colEnd = &matColIdx_h(matRowPtr_h(irow + 1));

      for (size_t jn = 0; jn < nodeList.size(); ++jn) {
        auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
        matValues_h(matRowPtr_h(irow) + pos) += kele[in + jn * nodeList.size()];
      }
    }
  }

  printf("Processed %d MFEM_L elements.\n", mfemCount);
}

}  // namespace IMSI
