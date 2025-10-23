#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <functional>
#include <optional>

#include "Element.h"
#include "MathUtils.h"
#include "MeshUtils.h"
#include "QuadratureRule.h"
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

  template <typename Device>
  void
  GetLinearSystem_v(
      Kokkos::View<double*, Device> rhs,
      Kokkos::View<size_t*, Device> matRowPtr,
      Kokkos::View<int*, Device>    matColIdx,
      Kokkos::View<double*, Device> matValues) const;

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

  template <int dim, int nNodes, typename ElementClass>
  void
  ElementaryDataLagrangeFE_v(
      ElementClass&                              element,
      const std::array<simd_type, dim * nNodes>& coords_v,
      simd_type*                                 rele,
      simd_type*                                 kele) const
  {
    //
    std::array<simd_type, dim*(dim + 1)> pointJac;
    std::array<simd_type, dim>           alpha;
    std::array<simd_type, nNodes * dim>  GradPhi;
    //
    for (int iq = 0; iq < ruleLength; ++iq) {
      auto NandGradN = element.GetValuesGradients(simd_type(xi[iq]), simd_type(eta[iq]), simd_type(zeta[iq]));
      pointJac.fill(simd_type(0));
      for (int jd = 0; jd <= dim; ++jd) {
        for (int id = 0; id < dim; ++id) {
          for (int kn = 0; kn < nNodes; ++kn) {
            pointJac[id + jd * dim] += NandGradN[kn + jd * nNodes] * coords_v[id + kn * dim];
          }
        }
      }
      auto const xq = pointJac[0];
      auto const yq = (dim > 1) ? pointJac[1] : simd_type(0);
      auto const zq = (dim > 2) ? pointJac[2] : simd_type(0);
      if (ax.has_value()) { alpha[0] = SIMDize(xq, yq, zq, ax.value()); }
      if constexpr (dim > 1) {
        if (ay.has_value()) { alpha[1] = SIMDize(xq, yq, zq, ay.value()); }
      }
      if constexpr (dim > 2) {
        if (az.has_value()) { alpha[2] = SIMDize(xq, yq, zq, az.value()); }
      }
      //
      // Get the inverse of the Jacobian
      //
      simd_type detJ(1);
      auto* __restrict J = &pointJac[dim];
      InverseInPlace<dim>(J, detJ);
      //
      auto* __restrict GradN = &NandGradN[nNodes];
      GradPhi.fill(simd_type(0));
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in < dim; ++in) {
          for (int kn = 0; kn < dim; ++kn) { GradPhi[in + jn * dim] += J[in + kn * dim] * GradN[jn + kn * nNodes]; }
        }
      }
      //
      simd_type w_v(weight[iq]);
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = 0; in <= jn; ++in) {
          for (int kn = 0; kn < dim; ++kn) {
            kele[in + jn * nNodes] += GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim] * w_v * detJ;
          }
        }
      }
      // Symmetrize the matrix
      for (int jn = 0; jn < nNodes; ++jn) {
        for (int in = jn + 1; in < nNodes; ++in) { kele[in + jn * nNodes] = kele[jn + in * nNodes]; }
      }
      //
      simd_type fq(0);
      if (f.has_value()) { fq = SIMDize(xq, yq, zq, f.value()); }
      for (int in = 0; in < nNodes; ++in) { rele[in] += fq * NandGradN[in] * w_v * detJ; }
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
        Kokkos::TeamPolicy<Device>(Kokkos::num_threads(), 1),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type& team) {
          auto const          ratio = (eleList.length + team.league_size() - 1) / team.league_size();
          auto const          kptr  = team.team_rank() + ratio * team.league_rank();
          std::vector<double> rele(maxNumDofsPerEle);
          std::vector<double> kele(maxNumDofsPerEle * maxNumDofsPerEle);
          for (int jj = 0; jj < ratio; ++jj) {
            auto const ik = kptr + jj;
            if (ik < eleList.length) {
              auto const eleID    = eleList(ik);
              auto       nodeList = meshInfo.mesh.NodeList(eleID);
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
                      fe1DQ1 element;
                      this->ElementaryDataLagrangeFE<1, 2, fe1DQ1>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                    default:
                    case 2: {
                      fe2DQ1 element;
                      this->ElementaryDataLagrangeFE<2, 4, fe2DQ1>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                    case 3: {
                      fe3DQ1 element;
                      this->ElementaryDataLagrangeFE<3, 8, fe3DQ1>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                  }
                  break;
                }
                case ElementType::Q2: {
                  switch (sdim) {
                    case 1: {
                      fe1DQ2 element;
                      this->ElementaryDataLagrangeFE<1, 3, fe1DQ2>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                    default:
                    case 2: {
                      fe2DQ2 element;
                      this->ElementaryDataLagrangeFE<2, 9, fe2DQ2>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                    case 3: {
                      fe3DQ2 element;
                      this->ElementaryDataLagrangeFE<3, 27, fe3DQ2>(element, nodeList, &rele[0], &kele[0]);
                      break;
                    }
                  }
                }
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
            }
          }
        });
    Kokkos::fence();
  }
}

template <typename Device>
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

  constexpr int vecSize = Kokkos::Experimental::native_simd<double>::size();

  for (int ic = 0; ic < c2e.numRows(); ++ic) {
    auto const eleList = c2e.rowConst(ic);
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<Device>(Kokkos::num_threads(), 1),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type& team) {
          auto const             klen = (eleList.length + team.league_size() - 1) / team.league_size();
          auto const             kptr = team.team_rank() + klen * team.league_rank();
          std::vector<simd_type> rele_v(maxNumDofsPerEle);
          std::vector<simd_type> kele_v(maxNumDofsPerEle * maxNumDofsPerEle);
          int                    jj = 0;
          for (; jj + vecSize <= klen; jj += vecSize) {
            //
            // Element type for eleID
            // !!! WARNING !!! This step assumes that the 'vecSize' elements are of the same type.
            //
            auto const ik       = kptr + jj;
            if (ik + vecSize > eleList.length) {
              break;
            }
            auto const eleID    = eleList(ik);
            int        numNodes = 0;
            switch (meshInfo.mesh.GetCellType(eleID)) {
              default:
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
                    rele_v.assign(numNodes, simd_type(0));
                    kele_v.assign(numNodes * numNodes, simd_type(0));
                    //
                    std::array<simd_type, fe2DQ1::numNode * fe2DQ1::sdim> coords_v{};
                    {
                      std::array<double, fe2DQ1::numNode * fe2DQ1::sdim * vecSize> coords{};
                      for (int jE = 0; jE < vecSize; ++jE) {
                        auto const nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                        for (int i = 0; i < fe2DQ1::numNode; ++i) {
                          auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                          {
                              auto const shift = jE + i * fe2DQ1::sdim * vecSize;
                            coords[shift] = vertex[0];
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
                    this->ElementaryDataLagrangeFE_v<fe2DQ1::sdim, fe2DQ1::numNode, fe2DQ1>(
                        element, coords_v, &rele_v[0], &kele_v[0]);
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
                    rele_v.assign(numNodes, simd_type(0));
                    kele_v.assign(numNodes * numNodes, simd_type(0));
                    //
                    std::array<simd_type, fe2DQ2::numNode * fe2DQ2::sdim> coords_v{};
                    {
                      std::array<double, coords_v.size() * vecSize> coords{};
                      for (int jE = 0; jE < vecSize; ++jE) {
                        auto const nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                        for (int i = 0; i < fe2DQ2::numNode; ++i) {
                          auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                          {
                            auto const shift = jE + i * fe2DQ2::sdim * vecSize;
                            coords[shift] = vertex[0];
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
                    this->ElementaryDataLagrangeFE_v<fe2DQ2::sdim, fe2DQ2::numNode, fe2DQ2>(
                        element, coords_v, &rele_v[0], &kele_v[0]);
                    break;
                  }
                  case 3: {
                    // fe3DQ2 element;
                    // this->ElementaryDataLagrangeFE<3, 27, fe3DQ2>(element, nodeList, &rele[0],
                    //                                               &kele[0]);
                    // break;
                  }
                }
              }
            }
            //
            {
              std::vector<double> rele(numNodes * vecSize);
              for (int i = 0; i < numNodes; ++i) {
                rele_v[i].copy_to(&rele[i * vecSize], Kokkos::Experimental::element_aligned_tag());
              }
              for (int jE = 0; jE < vecSize; ++jE) {
                auto const& nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                for (int in = 0; in < size(nodeList); ++in) { rhs(nodeList[in]) += rele[jE + in * vecSize]; }
              }
            }
            //
            {
              std::vector<double> kele(numNodes * numNodes * vecSize);
              for (int i = 0; i < numNodes * numNodes; ++i) {
                kele_v[i].copy_to(&kele[i * vecSize], Kokkos::Experimental::element_aligned_tag());
              }
              for (int jE = 0; jE < vecSize; ++jE) {
                auto const& nodeList = meshInfo.mesh.NodeList(eleList(ik + jE));
                for (int in = 0; in < size(nodeList); ++in) {
                  auto const irow     = nodeList[in];
                  auto const colBegin = &matColIdx(matRowPtr(irow));
                  auto const colEnd   = &matColIdx(matRowPtr(irow + 1));
                  for (int jn = 0; jn < size(nodeList); ++jn) {
                    auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
                    matValues(matRowPtr(irow) + pos) += kele[jE + in * vecSize + jn * numNodes * vecSize];
                  }
                }
              }
            }
          }
          //
          // Treat the remaining elements if needed
          //
          std::vector<double> rele(maxNumDofsPerEle);
          std::vector<double> kele(maxNumDofsPerEle * maxNumDofsPerEle);
          for (; jj < klen; jj += 1) {
            auto const ik       = kptr + jj;
            if (ik >= eleList.length) {
              break;
            }
            auto const eleID    = eleList(ik);
            auto       nodeList = meshInfo.mesh.NodeList(eleID);
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
                    fe1DQ1 element;
                    this->ElementaryDataLagrangeFE<1, 2, fe1DQ1>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                  default:
                  case 2: {
                    fe2DQ1 element;
                    this->ElementaryDataLagrangeFE<2, 4, fe2DQ1>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                  case 3: {
                    fe3DQ1 element;
                    this->ElementaryDataLagrangeFE<3, 8, fe3DQ1>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                }
                break;
              }
              case ElementType::Q2: {
                switch (sdim) {
                  case 1: {
                    fe1DQ2 element;
                    this->ElementaryDataLagrangeFE<1, 3, fe1DQ2>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                  default:
                  case 2: {
                    fe2DQ2 element;
                    this->ElementaryDataLagrangeFE<2, 9, fe2DQ2>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                  case 3: {
                    fe3DQ2 element;
                    this->ElementaryDataLagrangeFE<3, 27, fe3DQ2>(element, nodeList, &rele[0], &kele[0]);
                    break;
                  }
                }
              }
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
          }  // for (; jj < kLen; jj += 1)
        });
    Kokkos::fence();
  }  // for (int ic = 0; ic < ; ++ic)
}

}  // namespace IMSI
