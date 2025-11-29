#pragma once

#include <functional>
#include <optional>

#include "../main_config.h"
#include "Element.h"
#include "MathUtils.h"
#include "MeshUtils.h"
#include "QuadratureRule.h"
#include "SymmetricSparse.hpp"
#include "fe1DQ1.h"
#include "fe1DQ2.h"
#include "fe2DQ1.h"
#include "fe2DQ2.h"
#include "fe3DQ1.h"
#include "fe3DQ2.h"

namespace IMSI {

class ScaledLaplacian
{
 public:
  ScaledLaplacian(
      const MeshConnectivity<>&                     meshData,
      std::function<double(double, double, double)> alpha_x,
      std::function<double(double, double, double)> beta_y,
      std::function<double(double, double, double)> f_rhs,
      RuleType                                      quadRule,
      int                                           quadOrder)
      : meshInfo(meshData),
        ax(alpha_x),
        ay(beta_y),
        az(std::nullopt),
        f(f_rhs),
        ruleType(quadRule),
        ruleOrder(quadOrder)
  {
    auto const sdim = meshInfo.mesh.GetSpatialDimension();
    getQuadrature(ruleType, sdim, ruleOrder, ruleLength, weight, xi, eta, zeta);
  }

  template <typename Device>
  void
  GetLinearSystem(
      Kokkos::View<double*, Device> rhs,
      Kokkos::View<size_t*, Device> matRowPtr,
      Kokkos::View<int*, Device>    matColIdx,
      Kokkos::View<double*, Device> matValues) const;

  template <typename Device, bool useSIMD, bool useColoring = true>
  void
  GetLinearSystem_v(
      Kokkos::View<double*, Device> rhs,
      Kokkos::View<size_t*, Device> matRowPtr,
      Kokkos::View<int*, Device>    matColIdx,
      Kokkos::View<double*, Device> matValues) const;

mutable std::vector< std::vector<double> > phiMFEM;
  int const ratio            = 32;

  void OutputMFEMFine(double* uCoarse, int numEleX, int numEleY) const {
    //
    int numFineEleX = ratio * numEleX;
    int numFineEleY = ratio * numEleY;
    int numFineNodes = (numFineEleX + 1) * (numFineEleY + 1);
    std::vector<double> uFine(numFineNodes, 0);
    std::vector<double> uTrace(4);
    //
    for (int iy = 0; iy < numEleY; ++iy) {
      for (int ix = 0; ix < numEleX; ++ix) {
        int eleID = ix + iy * numEleX;
        auto const& phi = phiMFEM[eleID];
        uTrace[0] = uCoarse[ix + iy * (numEleX + 1)];
        uTrace[1] = uCoarse[ix + 1 + iy * (numEleX + 1)];
        uTrace[2] = uCoarse[ix + 1 + (iy + 1) * (numEleX + 1)];
        uTrace[3] = uCoarse[ix + (iy + 1) * (numEleX + 1)];
        for (int jy = 0; jy <= ratio; ++jy) {
          for (int jx = 0; jx <= ratio; ++jx) {
            int nodeID = ix * ratio + jx + (iy * ratio + jy) * (numFineEleX + 1);
            uFine[nodeID] = 0;
            for (int k = 0; k < 4; ++k) {
              uFine[nodeID] += phi[jx + jy * (ratio + 1) + k * (ratio + 1) * (ratio + 1)] * uTrace[k];
            }
          }
        }
      }
    }
    //
    std::ofstream outFine("outputFine.txt");
    for (int i = 0; i < numFineNodes; ++i) {
      outFine << uFine[i] << std::endl;
    }
    outFine.close();
  }

 protected:
  const MeshConnectivity<>                                     meshInfo;
  std::optional<std::function<double(double, double, double)>> ax;
  std::optional<std::function<double(double, double, double)>> ay;
  std::optional<std::function<double(double, double, double)>> az;
  std::optional<std::function<double(double, double, double)>> f;

  RuleType ruleType   = RuleType::Gauss;
  int      ruleOrder  = 1;
  int      ruleLength = 0;

  std::vector<double> weight;
  std::vector<double> xi, eta, zeta;

 protected:
  using simd_type = Kokkos::Experimental::native_simd<double>;

  static simd_type
  SIMDize(simd_type x, simd_type y, simd_type z, const std::function<double(double, double, double)>& g)
  {
    simd_type                             val(0);
    std::array<double, simd_type::size()> xa, ya, za, ga;
    x.copy_to(&xa[0], Kokkos::Experimental::element_aligned_tag());
    y.copy_to(&ya[0], Kokkos::Experimental::element_aligned_tag());
    z.copy_to(&za[0], Kokkos::Experimental::element_aligned_tag());
    for (int i = 0; i < simd_type::size(); ++i) { ga[i] = g(xa[i], ya[i], za[i]); }
    val.copy_from(&ga[0], Kokkos::Experimental::element_aligned_tag());
    return val;
  }

  template <int dim, int nNodes, typename ElementClass, typename Scalar>
  void
  ElementaryDataLagrangeFE_impl(
      ElementClass&                           element,
      const std::array<Scalar, nNodes * dim>& nodes,
      Scalar                                  w,
      Scalar                                  xi,
      Scalar                                  eta,
      Scalar                                  zeta,
      Scalar*                                 rele,
      Scalar*                                 kele) const
  {
    std::array<Scalar, dim*(dim + 1)> pointJac;
    std::array<Scalar, dim>           alpha;
    std::array<Scalar, nNodes * dim>  GradPhi;
    //
    auto NandGradN = element.GetValuesGradients(xi, eta, zeta);
    pointJac.fill(Scalar(0));
    for (int jd = 0; jd <= dim; ++jd) {
      for (int id = 0; id < dim; ++id) {
        for (int kn = 0; kn < nNodes; ++kn) {
          pointJac[id + jd * dim] += NandGradN[kn + jd * nNodes] * nodes[id + kn * dim];
        }
      }
    }
    auto const xq = pointJac[0];
    auto const yq = (dim > 1) ? pointJac[1] : Scalar(0);
    auto const zq = (dim > 2) ? pointJac[2] : Scalar(0);
    if (ax.has_value()) { alpha[0] = ax->operator()(xq, yq, zq); }
    if constexpr (dim > 1) {
      if (ay.has_value()) { alpha[1] = ay->operator()(xq, yq, zq); }
    }
    if constexpr (dim > 2) {
      if (az.has_value()) { alpha[2] = az->operator()(xq, yq, zq); }
    }
    //
    // Get the inverse of the Jacobian
    //
    Scalar detJ(1);
    Scalar* __restrict J = &pointJac[dim];
    InverseInPlace<dim>(J, detJ);
    //
    Scalar* __restrict GradN = &NandGradN[nNodes];
    GradPhi.fill(Scalar(0));
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = 0; in < dim; ++in) {
        for (int kn = 0; kn < dim; ++kn) { GradPhi[in + jn * dim] += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
      }
    }
    //
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = 0; in <= jn; ++in) {
        for (int kn = 0; kn < dim; ++kn) {
          kele[in + jn * nNodes] += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim] * w * detJ;
        }
      }
    }
    // Symmetrize the matrix
    for (int jn = 0; jn < nNodes; ++jn) {
      for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
    }
    //
    Scalar fq(0);
    if (f.has_value()) { fq = f->operator()(xq, yq, zq); }
    for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * w * detJ; }
  }

  template <int dim, int nNodes, typename ElementClass>
  void
  ElementaryDataLagrangeFE(ElementClass& element, const std::vector<int>& nodeList, double* rele, double* kele) const
  {
    std::array<double, nNodes * dim> nodes;
    for (int i = 0; i < nNodes; ++i) {
      auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
      std::copy(&vertex[0], &vertex[0] + dim, &nodes[i * dim]);
    }
    //
    std::array<double, dim*(dim + 1)> pointJac;
    std::array<double, dim>           alpha;
    std::array<double, nNodes * dim>  GradPhi;
    //
    for (int iq = 0; iq < ruleLength; ++iq) {
      auto NandGradN = element.GetValuesGradients(xi[iq], eta[iq], zeta[iq]);
      pointJac.fill(0.0);
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          for (int kn = 0; kn < nNodes; ++kn) {
            pointJac[id + jd * dim] += NandGradN[kn + jd * nNodes] * nodes[id + kn * dim];
          }
        }
      }
      auto const xq = pointJac[0];
      auto const yq = (dim > 1) ? pointJac[1] : double(0.0);
      auto const zq = (dim > 2) ? pointJac[2] : double(0.0);
      if (ax.has_value()) { alpha[0] = ax->operator()(xq, yq, zq); }
      if constexpr (dim > 1) {
        if (ay.has_value()) { alpha[1] = ay->operator()(xq, yq, zq); }
      }
      if constexpr (dim > 2) {
        if (az.has_value()) { alpha[2] = az->operator()(xq, yq, zq); }
      }
      //
      // Get the inverse of the Jacobian
      //
      double detJ          = 1.0;
      double* __restrict J = &pointJac[dim];
      InverseInPlace<dim>(J, detJ);
      //
      double* __restrict GradN = &NandGradN[nNodes];
      GradPhi.fill(0);
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          for (int kn = 0; kn < dim; ++kn) { GradPhi[in + jn * dim] += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
        }
      }
      //
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          for (int kn = 0; kn < dim; ++kn) {
            kele[in + jn * nNodes] += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim] * weight[iq] * detJ;
          }
        }
      }
      // Symmetrize the matrix
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
      }
      //
      double fq = 0.0;
      if (f.has_value()) { fq = f->operator()(xq, yq, zq); }
      for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * weight[iq] * detJ; }
    }
  }

  template <typename Scalar, typename ElementClass>
  void
  ElementaryDataLagrangeFE_t(
      ElementClass& element,
      const Scalar* __restrict__ coords_v,
      Scalar* __restrict__ rele,
      Scalar* __restrict__ kele) const
  {
    auto constexpr dim    = ElementClass::sdim;
    auto constexpr nNodes = ElementClass::numNode;
    //
    std::array<Scalar, dim*(dim + 1)> pointJac;
    std::array<Scalar, dim>           alpha;
    std::array<Scalar, nNodes * dim>  GradPhi;
    //
    bool has_ax = ax.has_value();
    bool has_ay = ay.has_value();
    bool has_az = az.has_value();
    bool has_f  = f.has_value();
    //
    for (int iq = 0; iq < ruleLength; ++iq) {
      auto NandGradN = element.GetValuesGradients(Scalar(xi[iq]), Scalar(eta[iq]), Scalar(zeta[iq]));
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          Scalar jacEntry{0};
          for (int kn = 0; kn < nNodes; ++kn) { jacEntry += NandGradN[kn + jd * nNodes] * coords_v[id + kn * dim]; }
          pointJac[id + jd * dim] = jacEntry;
        }
      }
      auto const xq = pointJac[0];
      auto const yq = (dim > 1) ? pointJac[1] : Scalar(0);
      auto const zq = (dim > 2) ? pointJac[2] : Scalar(0);
      if (has_ax) {
        if constexpr (std::is_same_v<simd_type, Scalar>) {
          alpha[0] = SIMDize(xq, yq, zq, ax.value());
        } else {
          alpha[0] = ax->operator()(xq, yq, zq);
        }
      }
      if constexpr (dim > 1) {
        if (has_ay) {
          if constexpr (std::is_same_v<simd_type, Scalar>) {
            alpha[1] = SIMDize(xq, yq, zq, ay.value());
          } else {
            alpha[1] = ay->operator()(xq, yq, zq);
          }
        }
      }
      if constexpr (dim > 2) {
        if (has_az) {
          if constexpr (std::is_same_v<simd_type, Scalar>) {
            alpha[2] = SIMDize(xq, yq, zq, az.value());
          } else {
            alpha[2] = az->operator()(xq, yq, zq);
          }
        }
      }
      //
      // Get the inverse of the Jacobian
      //
      Scalar detJ(1);
      auto* __restrict J = &pointJac[dim];
      InverseInPlace<dim>(J, detJ);
      //
      auto* __restrict GradN = &NandGradN[nNodes];
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          Scalar tmpGrad{0};
          for (int kn = 0; kn < dim; ++kn) { tmpGrad += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
          GradPhi[in + jn * dim] = tmpGrad;
        }
      }
      //
      Scalar w_v(weight[iq]);
      Scalar coeff = w_v * detJ;
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          Scalar sum(0);
          for (int kn = 0; kn < dim; ++kn) { sum += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim]; }
          kele[in + jn * nNodes] += sum * coeff;
        }
      }
      // Symmetrize the matrix
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
      }
      //
      Scalar fq(0);
      if (has_f) {
        if constexpr (std::is_same_v<simd_type, Scalar>) {
          fq = SIMDize(xq, yq, zq, f.value());
        } else {
          fq = f->operator()(xq, yq, zq);
        }
        for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * coeff; }
      }
    }
  }

  template <typename Scalar, typename ElementClass>
  void
  ElementaryDataMFEM_t(ElementClass& element, const Scalar* coords_v, Scalar* rele, Scalar* kele,
    std::vector<Scalar>& phi) const
  {
    /// TODO Generalize these parameters
    size_t    maxNumDofsPerEle = 4;
    int const numNodes         = (ratio + 1) * (ratio + 1);
    //
    std::vector<Scalar> rhs(numNodes, 0);
    std::vector<Scalar> rFineEle(maxNumDofsPerEle, 0);
    std::vector<Scalar> kFineEle(maxNumDofsPerEle * maxNumDofsPerEle, 0);
    //
    double const hx = (coords_v[2] - coords_v[0]) / double(ratio);
    double const hy = (coords_v[7] - coords_v[1]) / double(ratio);
    //
    // Fill sparsity of local problem
    //
    Kokkos::Timer    timer;
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
    std::vector<Scalar> matValues(matColIdx.size(), 0);
    //
    // Get numerical entries for local sparse matrix
    //
    std::array<Scalar, fe2DQ1::numNode * fe2DQ1::sdim> coords{};
    using simd_t                                                        = Kokkos::Experimental::native_simd<Scalar>;
    constexpr int                                                 width = simd_t::size();
    std::array<simd_t, fe2DQ1::numNode * fe2DQ1::sdim>            coords_simd;
    std::array<simd_t, fe2DQ1::numNode>                           rFineEle_simd;
    std::array<simd_t, fe2DQ1::numNode * fe2DQ1::numNode>         kFineEle_simd;
    std::array<Scalar, fe2DQ1::numNode * width>                   rFineEle_v;
    std::array<Scalar, fe2DQ1::numNode * fe2DQ1::numNode * width> kFineEle_v;

    timer.reset();
    for (int iy = 0; iy < ratio; ++iy) {
      int ix = 0;
      for (; ix <= ratio - width; ix += width) {
        for (int d = 0; d < 8; ++d) {
          std::array<Scalar, width> tmp_coords;
          for (int k = 0; k < width; ++k) {
            Scalar const cx = coords_v[0] + (ix + k) * hx;
            Scalar const cy = coords_v[1] + iy * hy;
            if (d == 0 || d == 6)
              tmp_coords[k] = cx;
            else if (d == 1 || d == 3)
              tmp_coords[k] = cy;
            else if (d == 2 || d == 4)
              tmp_coords[k] = cx + hx;
            else if (d == 5 || d == 7)
              tmp_coords[k] = cy + hy;
          }
          coords_simd[d].copy_from(tmp_coords.data(), Kokkos::Experimental::element_aligned_tag());
        }

        for (auto& val : rFineEle_simd) val = 0;
        for (auto& val : kFineEle_simd) val = 0;

        this->ElementaryDataLagrangeFE_t<simd_t>(element, &coords_simd[0], &rFineEle_simd[0], &kFineEle_simd[0]);

        for (int k = 0; k < fe2DQ1::numNode; ++k)
          rFineEle_simd[k].copy_to(&rFineEle_v[k * width], Kokkos::Experimental::element_aligned_tag());
        for (int k = 0; k < fe2DQ1::numNode * fe2DQ1::numNode; ++k)
          kFineEle_simd[k].copy_to(&kFineEle_v[k * width], Kokkos::Experimental::element_aligned_tag());

        for (int k = 0; k < width; ++k) {
          int const                        iix = ix + k;
          std::array<int, fe2DQ1::numNode> nodeList{
              iix + iy * (ratio + 1),
              iix + 1 + iy * (ratio + 1),
              iix + 1 + (iy + 1) * (ratio + 1),
              iix + (iy + 1) * (ratio + 1)};
          //
          for (int in = 0; in < size(nodeList); ++in) { rhs[nodeList[in]] += rFineEle_v[in * width + k]; }
          //
          for (int in = 0; in < size(nodeList); ++in) {
            auto const irow     = nodeList[in];
            auto const colBegin = &matColIdx[matRowPtr[irow]];
            auto const colEnd   = &matColIdx[matRowPtr[irow + 1]];
            for (int jn = 0; jn < size(nodeList); ++jn) {
              auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
              matValues[matRowPtr[irow] + pos] += kFineEle_v[(in + jn * size(nodeList)) * width + k];
            }
          }
        }
      }
      for (; ix < ratio; ++ix) {
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
        this->ElementaryDataLagrangeFE_t<Scalar>(element, &coords[0], &rFineEle[0], &kFineEle[0]);
        //
        for (int in = 0; in < size(nodeList); ++in) { rhs[nodeList[in]] += rFineEle[in]; }
        //
        for (int in = 0; in < size(nodeList); ++in) {
          auto const irow     = nodeList[in];
          auto const colBegin = &matColIdx[matRowPtr[irow]];
          auto const colEnd   = &matColIdx[matRowPtr[irow + 1]];
          for (int jn = 0; jn < size(nodeList); ++jn) {
            auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
            matValues[matRowPtr[irow] + pos] += kFineEle[in + jn * size(nodeList)];
          }
        }
      }
    }
    printf(" --- Local matrix assembly %e ", timer.seconds());
    //
    //--- Get map from global numbering to free numbering (and vice versa)
    //
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
    //
    int const           numVectors = 4;
    phi.resize(numVectors * numNodes, 0);
    for (int iy = 0; iy <= ratio; ++iy) {
      int    ix              = 0;
      int    in              = ix + iy * (ratio + 1);
      Scalar eta             = Scalar(iy) / Scalar(ratio);
      phi[in]                = Scalar(1) - eta;
      phi[in + 3 * numNodes] = eta;
      //
      ix                     = ratio;
      in                     = ix + iy * (ratio + 1);
      phi[in + numNodes]     = Scalar(1) - eta;
      phi[in + 2 * numNodes] = eta;
    }
    for (int ix = 0; ix <= ratio; ++ix) {
      int    iy          = 0;
      int    in          = ix + iy * (ratio + 1);
      Scalar xi          = Scalar(ix) / Scalar(ratio);
      phi[in]            = Scalar(1) - xi;
      phi[in + numNodes] = xi;
      //
      iy                     = ratio;
      in                     = ix + iy * (ratio + 1);
      phi[in + 3 * numNodes] = Scalar(1) - xi;
      phi[in + 2 * numNodes] = xi;
    }
    //
    //--- Prepare the linear system without boundary conditions
    //
    auto const       n = freeToGlobal.size();
    std::vector<int> newRowPtr(n + 1);
    newRowPtr[0] = 0;
    int count    = 0;
    for (int i = 0; i < n; ++i) {
      auto gDof   = freeToGlobal[i];
      int  iCount = 0;
      for (auto k = matRowPtr[gDof]; k < matRowPtr[gDof + 1]; ++k) { iCount += (globalToFree[matColIdx[k]] != -1); }
      newRowPtr[i + 1] = iCount;
      count += iCount;
    }
    //
    for (int i = 0; i < newRowPtr.size() - 1; ++i) { newRowPtr[i + 1] += newRowPtr[i]; }
    auto const newNNZ = newRowPtr[n];
    //
    std::vector<int>    newColIdx(newNNZ);
    std::vector<Scalar> newValues(newNNZ);
    // std::vector<Scalar> newRHS(n);
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
    //
    // Compute the basis functions
    //
    std::vector<Scalar> btmp(numVectors * n, 0);
    SparseMatrix<Scalar> K(
        numNodes, numNodes, matColIdx.size(), matRowPtr.data(), matColIdx.data(), matValues.data());
    std::vector<Scalar> Kphi(numVectors * numNodes);
    {
      K.Apply(numVectors, &(phi[0]), &Kphi[0]);
      for (int ii = 0; ii < n; ++ii) {
        for (int ir = 0; ir < numVectors; ++ir) { btmp[ii + ir * n] = -Kphi[freeToGlobal[ii] + ir * numNodes]; }
      }
    }
    //
    timer.reset();
    SymmetricSparse<Scalar> Ktmp(n, newNNZ, newRowPtr.data(), newColIdx.data(), newValues.data());
    Ktmp.factor();
    printf(" --- Local stiffness facto %e \n", timer.seconds());
    //
    std::vector<Scalar> utmp(numVectors * n, 0);
    Ktmp.Solve(numVectors, &btmp[0], &utmp[0]);
    for (int ii = 0; ii < n; ++ii) {
      int grow = freeToGlobal[ii];
      for (int ir = 0; ir < numVectors; ++ir) { phi[grow + ir * numNodes] = utmp[ii + ir * n]; }
    }
    //
    for (int ir = 0; ir < numVectors; ++ir) {
      Scalar sum = 0;
      for (int ii = 0; ii < numNodes; ++ii) { sum += phi[ii + ir * numNodes] * rhs[ii]; }
      rele[ir] = sum;
    }
    //
    {
      K.Apply(numVectors, &(phi[0]), &Kphi[0]);
      for (int ir = 0; ir < numVectors; ++ir) {
        for (int jr = 0; jr <= ir; ++jr) {
          Scalar sum = 0;
          for (int ii = 0; ii < numNodes; ++ii) { sum += phi[ii + ir * numNodes] * Kphi[ii + jr * numNodes]; }
          kele[ir + jr * numVectors] = sum;
          kele[jr + ir * numVectors] = sum;
        }
      }
    }
  }
};

}  // namespace IMSI

//
// Definition of functions
//

namespace IMSI {

template <typename Device>
void
ScaledLaplacian::GetLinearSystem(
    Kokkos::View<double*, Device> rhs,
    Kokkos::View<size_t*, Device> matRowPtr,
    Kokkos::View<int*, Device>    matColIdx,
    Kokkos::View<double*, Device> matValues) const
{
  size_t maxNumDofsPerEle = 0;
  Kokkos::parallel_reduce(
      "MaxDofsPerEle",
      Kokkos::RangePolicy<Device>(0, meshInfo.mesh.NumberCells()),
      KOKKOS_LAMBDA(const int& i, size_t& nMax) { nMax = std::max<size_t>(nMax, size(meshInfo.mesh.NodeList(i))); },
      Kokkos::Max<size_t>(maxNumDofsPerEle));

  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const eleList = c2e.rowConst(ic);
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<Device>(eleList.length, 1),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type& team) {
          auto const          ik = team.league_rank();
          std::vector<double> rele(maxNumDofsPerEle);
          std::vector<double> kele(maxNumDofsPerEle * maxNumDofsPerEle);
          auto const          eleID    = eleList(ik);
          auto                nodeList = meshInfo.mesh.NodeList(eleID);
          //
          rele.assign(size(nodeList), 0);
          kele.assign(size(nodeList) * size(nodeList), 0);
          //
          // Element type for eleID
          //
          switch (meshInfo.mesh.GetCellType(eleID)) {
            default:
            case ElementType::Q1: {
              switch (sdim) {
                case 1: {
                  std::array<double, fe1DQ1::numNode * fe1DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe1DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe1DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                default:
                case 2: {
                  std::array<double, fe2DQ1::numNode * fe2DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe2DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe2DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                case 3: {
                  std::array<double, fe3DQ1::numNode * fe3DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe3DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe3DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
              }
              break;
            }  // case ElementType::Q1:
            case ElementType::Q2: {
              switch (sdim) {
                case 1: {
                  std::array<double, fe1DQ2::numNode * fe1DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe1DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe1DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                default:
                case 2: {
                  std::array<double, fe2DQ2::numNode * fe2DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe2DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe2DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                case 3: {
                  std::array<double, fe3DQ2::numNode * fe3DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe3DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe3DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
              }
              break;
            }  // case ElementType::Q2
            case ElementType::MFEM_L: {
              switch (sdim) {
                default: {
                  exit(EXIT_FAILURE);
                  break;
                }
              }
              break;
            }  // case ElementType::MFEM_L:
          }
          //
          for (int in = 0; in < size(nodeList); ++in) { rhs(nodeList[in]) += rele[in]; }
          //
          for (int in = 0; in < size(nodeList); ++in) {
            auto const irow     = nodeList[in];
            auto const colBegin = &matColIdx(matRowPtr(irow));
            auto const colEnd   = &matColIdx(matRowPtr(irow + 1));
            for (int jn = 0; jn < size(nodeList); ++jn) {
              auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
              matValues(matRowPtr(irow) + pos) += kele[in + jn * size(nodeList)];
            }
          }
        });
    Kokkos::fence();
  }
}

template <typename Device, bool useSIMD, bool useColoring>
void
ScaledLaplacian::GetLinearSystem_v(
    Kokkos::View<double*, Device> rhs,
    Kokkos::View<size_t*, Device> matRowPtr,
    Kokkos::View<int*, Device>    matColIdx,
    Kokkos::View<double*, Device> matValues) const
{
  size_t maxNumDofsPerEle = 0;
  Kokkos::parallel_reduce(
      "MaxDofsPerEle",
      Kokkos::RangePolicy<Device>(0, meshInfo.mesh.NumberCells()),
      KOKKOS_LAMBDA(const int& i, size_t& nMax) { nMax = std::max<size_t>(nMax, size(meshInfo.mesh.NodeList(i))); },
      Kokkos::Max<size_t>(maxNumDofsPerEle));

  auto const& c2e  = meshInfo.c2e;
  auto const  sdim = meshInfo.mesh.GetSpatialDimension();

  // TODO Remove the mutable characteristic
  phiMFEM.resize(meshInfo.mesh.NumberCells());

  int constexpr vecSize = (useSIMD) ? Kokkos::Experimental::native_simd<double>::size() : 1;
  printf(" --- Vectorized Assembly --- %d \n", vecSize);
  int const numColors = (useColoring) ? c2e.numRows() : 1;

  // Calculate scratch memory size for scalar element assembly
  const size_t scratchSize_rele = maxNumDofsPerEle * sizeof(double);
  const size_t scratchSize_kele = maxNumDofsPerEle * maxNumDofsPerEle * sizeof(double);
  const size_t totalScratchSize = scratchSize_rele + scratchSize_kele;

  for (int ic = 0; ic < numColors; ++ic) {
    [[maybe_unused]] auto const eleList   = c2e.rowConst(ic);
    auto const                  numEle    = (useColoring) ? eleList.length : meshInfo.mesh.NumberCells();
    auto const                  offSetLen = (vecSize == 1) ? numEle : numEle % vecSize;
    Kokkos::parallel_for(
        "ScalarAssembly",
        Kokkos::TeamPolicy<Device>(offSetLen, Kokkos::AUTO)
          .set_scratch_size(0, Kokkos::PerThread(totalScratchSize)),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type& team) {
          auto const          ik = team.league_rank();

          // Use scratch memory instead of heap allocation
          double* rele = (double*)team.thread_scratch(0).get_shmem(scratchSize_rele);
          double* kele = (double*)team.thread_scratch(0).get_shmem(scratchSize_kele);

          auto const          eleID    = (useColoring) ? eleList(ik) : ik;
          auto                nodeList = meshInfo.mesh.NodeList(eleID);
          //
          // Initialize scratch memory to zero
          for (int i = 0; i < size(nodeList); ++i) { rele[i] = 0.0; }
          for (int i = 0; i < size(nodeList) * size(nodeList); ++i) { kele[i] = 0.0; }
          //
          // Element type for eleID
          //
          switch (meshInfo.mesh.GetCellType(eleID)) {
            default:
            case ElementType::Q1: {
              switch (sdim) {
                case 1: {
                  std::array<double, fe1DQ1::numNode * fe1DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe1DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe1DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                default:
                case 2: {
                  std::array<double, fe2DQ1::numNode * fe2DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe2DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe2DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                case 3: {
                  std::array<double, fe3DQ1::numNode * fe3DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe3DQ1 element;
                  this->ElementaryDataLagrangeFE_t<double, fe3DQ1>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
              }
              break;
            }  // case ElementType::Q1:
            case ElementType::Q2: {
              switch (sdim) {
                case 1: {
                  std::array<double, fe1DQ2::numNode * fe1DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe1DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe1DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                default:
                case 2: {
                  std::array<double, fe2DQ2::numNode * fe2DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe2DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe2DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
                case 3: {
                  std::array<double, fe3DQ2::numNode * fe3DQ2::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe3DQ2 element;
                  this->ElementaryDataLagrangeFE_t<double, fe3DQ2>(element, &coords[0], &rele[0], &kele[0]);
                  break;
                }
              }
              break;
            }  // case ElementType::Q2
            case ElementType::MFEM_L: {
              switch (sdim) {
                default: {
                  exit(EXIT_FAILURE);
                  break;
                }
                case 2: {
                  std::array<double, fe2DQ1::numNode * fe2DQ1::sdim> coords{};
                  for (int i = 0; i < size(nodeList); ++i) {
                    auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                    std::copy(&vertex[0], &vertex[0] + sdim, &coords[i * sdim]);
                  }
                  fe2DQ1 element;
                  this->ElementaryDataMFEM_t<double>(element, &coords[0], &rele[0], &kele[0], phiMFEM[eleID]);
                  break;
                }
              }
              break;
            }  // case ElementType::MFEM_L:
          }
          //
          for (int in = 0; in < size(nodeList); ++in) {
            if (useColoring) {
              rhs(nodeList[in]) += rele[in];
            } else {
              Kokkos::atomic_add(&rhs(nodeList[in]), rele[in]);
            }
          }
          //
          for (int in = 0; in < size(nodeList); ++in) {
            auto const irow     = nodeList[in];
            auto const colBegin = &matColIdx(matRowPtr(irow));
            auto const colEnd   = &matColIdx(matRowPtr(irow + 1));
            for (int jn = 0; jn < size(nodeList); ++jn) {
              auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
              if (useColoring) {
                matValues(matRowPtr(irow) + pos) += kele[in + jn * size(nodeList)];
              } else {
                Kokkos::atomic_add(&matValues(matRowPtr(irow) + pos), kele[in + jn * size(nodeList)]);
              }
            }
          }
        });
    //
    // Calculate scratch memory size for SIMD element assembly
    const size_t scratchSize_rele_v = maxNumDofsPerEle * sizeof(simd_type);
    const size_t scratchSize_kele_v = maxNumDofsPerEle * maxNumDofsPerEle * sizeof(simd_type);
    const size_t scratchSize_coords = sdim * vecSize * maxNumDofsPerEle * sizeof(double);
    const size_t scratchSize_rele_unpacked = maxNumDofsPerEle * vecSize * sizeof(double);
    const size_t scratchSize_kele_unpacked = maxNumDofsPerEle * maxNumDofsPerEle * vecSize * sizeof(double);
    const size_t totalScratchSize_SIMD = scratchSize_rele_v + scratchSize_kele_v + scratchSize_coords +
                                          scratchSize_rele_unpacked + scratchSize_kele_unpacked;

    Kokkos::parallel_for(
        "SIMDAssembly",
        Kokkos::TeamPolicy<Device>((numEle - offSetLen) / vecSize, Kokkos::AUTO)
          .set_scratch_size(0, Kokkos::PerThread(totalScratchSize_SIMD)),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type& team) {
          // Use scratch memory instead of heap allocation
          simd_type* rele_v = (simd_type*)team.thread_scratch(0).get_shmem(scratchSize_rele_v);
          simd_type* kele_v = (simd_type*)team.thread_scratch(0).get_shmem(scratchSize_kele_v);
          double* coords = (double*)team.thread_scratch(0).get_shmem(scratchSize_coords);
          double* rele_unpacked = (double*)team.thread_scratch(0).get_shmem(scratchSize_rele_unpacked);
          double* kele_unpacked = (double*)team.thread_scratch(0).get_shmem(scratchSize_kele_unpacked);
          //
          // Element type for eleID
          // !!! WARNING !!! This step assumes that the 'vecSize' elements are of the same type.
          //
          auto const ik       = offSetLen + team.league_rank() * vecSize;
          auto const eleID    = (useColoring) ? eleList(ik) : ik;
          int        numNodes = 0;
          switch (meshInfo.mesh.GetCellType(eleID)) {
            default: {
              exit(EXIT_FAILURE);
              break;
            }
            case ElementType::Q1: {
              switch (sdim) {
                case 1: {
                  // fe1DQ1 element;
                  // this->ElementaryDataLagrangeFE<1, 2, fe1DQ1>(element, nodeList, &rele[0],
                  //                                              &kele[0]);
                  // break;
                }
                default:
                case 2: {
                  //
                  numNodes = fe2DQ1::numNode;
                  // Initialize scratch memory to zero
                  for (int i = 0; i < numNodes; ++i) { rele_v[i] = simd_type(0); }
                  for (int i = 0; i < numNodes * numNodes; ++i) { kele_v[i] = simd_type(0); }
                  //
                  std::array<simd_type, fe2DQ1::numNode * fe2DQ1::sdim> coords_v{};
                  {
                    for (int jE = 0; jE < vecSize; ++jE) {
                      auto const nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                      for (int i = 0; i < fe2DQ1::numNode; ++i) {
                        auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                        {
                          auto const shift        = jE + i * fe2DQ1::sdim * vecSize;
                          coords[shift]           = vertex[0];
                          coords[shift + vecSize] = vertex[1];
                        }
                      }
                    }
                    for (int i = 0; i < fe2DQ1::numNode * fe2DQ1::sdim; ++i) {
                      coords_v[i].copy_from(&coords[i * vecSize], Kokkos::Experimental::element_aligned_tag());
                    }
                  }
                  //
                  fe2DQ1 element;
                  this->ElementaryDataLagrangeFE_t<simd_type, fe2DQ1>(element, &coords_v[0], &rele_v[0], &kele_v[0]);
                  break;
                }
                case 3: {
                  // fe3DQ1 element;
                  // this->ElementaryDataLagrangeFE<3, 8, fe3DQ1>(element, nodeList, &rele[0],
                  //                                              &kele[0]);
                  // break;
                }
              }
              break;
            }
            case ElementType::Q2: {
              switch (sdim) {
                case 1: {
                  // fe1DQ2 element;
                  // this->ElementaryDataLagrangeFE<1, 3, fe1DQ2>(element, nodeList, &rele[0],
                  //                                              &kele[0]);
                  // break;
                }
                default:
                case 2: {
                  numNodes = fe2DQ2::numNode;
                  // Initialize scratch memory to zero
                  for (int i = 0; i < numNodes; ++i) { rele_v[i] = simd_type(0); }
                  for (int i = 0; i < numNodes * numNodes; ++i) { kele_v[i] = simd_type(0); }
                  //
                  std::array<simd_type, fe2DQ2::numNode * fe2DQ2::sdim> coords_v{};
                  {
                    for (int jE = 0; jE < vecSize; ++jE) {
                      auto const nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                      for (int i = 0; i < fe2DQ2::numNode; ++i) {
                        auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                        {
                          auto const shift        = jE + i * fe2DQ2::sdim * vecSize;
                          coords[shift]           = vertex[0];
                          coords[shift + vecSize] = vertex[1];
                        }
                      }
                    }
                    for (int i = 0; i < coords_v.size(); ++i) {
                      coords_v[i].copy_from(&coords[i * vecSize], Kokkos::Experimental::element_aligned_tag());
                    }
                  }
                  //
                  fe2DQ2 element;
                  this->ElementaryDataLagrangeFE_t<simd_type, fe2DQ2>(element, &coords_v[0], &rele_v[0], &kele_v[0]);
                  break;
                }
                case 3: {
                  // fe3DQ2 element;
                  //// TO BE IMPLEMENTED
                  break;
                }
              }
              break;
            }  // case ElementType::Q2:
          }
          //
          {
            // Use scratch memory for unpacking (already allocated above)
            for (int i = 0; i < numNodes; ++i) {
              rele_v[i].copy_to(&rele_unpacked[i * vecSize], Kokkos::Experimental::element_aligned_tag());
            }
            for (int jE = 0; jE < vecSize; ++jE) {
              auto const& nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
              for (int in = 0; in < size(nodeList); ++in) {
                if (useColoring) {
                  rhs(nodeList[in]) += rele_unpacked[jE + in * vecSize];
                } else {
                  Kokkos::atomic_add(&rhs(nodeList[in]), rele_unpacked[jE + in * vecSize]);
                }
              }
            }
          }
          //
          {
            // Use scratch memory for unpacking (already allocated above)
            for (int i = 0; i < numNodes * numNodes; ++i) {
              kele_v[i].copy_to(&kele_unpacked[i * vecSize], Kokkos::Experimental::element_aligned_tag());
            }
            for (int jE = 0; jE < vecSize; ++jE) {
              auto const& nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
              for (int in = 0; in < size(nodeList); ++in) {
                auto const irow     = nodeList[in];
                auto const colBegin = &matColIdx(matRowPtr(irow));
                auto const colEnd   = &matColIdx(matRowPtr(irow + 1));
                for (int jn = 0; jn < size(nodeList); ++jn) {
                  auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
                  if (useColoring) {
                    matValues(matRowPtr(irow) + pos) += kele_unpacked[jE + in * vecSize + jn * numNodes * vecSize];
                  } else {
                    Kokkos::atomic_add(
                        &matValues(matRowPtr(irow) + pos), kele_unpacked[jE + in * vecSize + jn * numNodes * vecSize]);
                  }
                }
              }
            }
          }
        });
    Kokkos::fence();
  }  // for (int ic = 0; ic < ; ++ic)
}

}  // namespace IMSI
