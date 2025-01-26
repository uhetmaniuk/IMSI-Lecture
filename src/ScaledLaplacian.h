#pragma once

#include "Element.h"
#include "fe1DQ1.h"
#include "fe1DQ2.h"
#include "fe2DQ1.h"
#include "fe2DQ2.h"
#include "fe3DQ1.h"
#include "fe3DQ2.h"
#include "MathUtils.h"
#include "MeshUtils.h"
#include "QuadratureRule.h"

#include <Kokkos_Core.hpp>

#include <functional>
#include <optional>

namespace IMSI {

    class ScaledLaplacian {
    public:
        ScaledLaplacian(const MeshConnectivity<> &meshData,
                        std::function<double(double, double, double)> alpha_x,
                        std::function<double(double, double, double)> beta_y,
                        std::function<double(double, double, double)> f_rhs,
                        RuleType quadRule, int quadOrder)
                : meshInfo(meshData), ax(alpha_x), ay(beta_y), az(std::nullopt), f(f_rhs),
                  ruleType(quadRule), ruleOrder(quadOrder) {
            auto const sdim = meshInfo.mesh.GetSpatialDimension();
            getQuadrature(ruleType, sdim, ruleOrder, ruleLength, weight, xi, eta, zeta);
        }

        template<typename Device>
        void GetLinearSystem(Kokkos::View<double *, Device> rhs,
                             Kokkos::View<size_t *, Device> matRowPtr,
                             Kokkos::View<int *, Device> matColIdx,
                             Kokkos::View<double *, Device> matValues) const;


    protected:

        const MeshConnectivity<> meshInfo;
        std::optional<std::function<double(double, double, double)> > ax;
        std::optional<std::function<double(double, double, double)> > ay;
        std::optional<std::function<double(double, double, double)> > az;
        std::optional<std::function<double(double, double, double)> > f;

        RuleType ruleType = RuleType::Gauss;
        int ruleOrder = 1;
        int ruleLength = 0;

        std::vector<double> weight;
        std::vector<double> xi, eta, zeta;

    protected:

        template<int dim, int nNodes, typename ElementClass>
        void ElementaryDataLagrangeFE(ElementClass &element,
                                      const std::vector<int> &nodeList,
                                      double *rele, double *kele) const {
            std::array<double, nNodes * dim> nodes;
            for (int i = 0; i < nNodes; ++i) {
                auto const vertex = meshInfo.mesh.GetVertex(nodeList[i]);
                std::copy(&vertex[0], &vertex[0] + dim, &nodes[i * dim]);
            }
            //
            std::array<double, dim * (dim + 1)> pointJac;
            std::array<double, dim> alpha;
            std::array<double, nNodes * dim> GradPhi;
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
                if (ax.has_value()) {
                    alpha[0] = ax->operator()(xq, yq, zq);
                }
                if constexpr (dim > 1) {
                    if (ay.has_value()) {
                        alpha[1] = ay->operator()(xq, yq, zq);
                    }
                }
                if constexpr (dim > 2) {
                    if (az.has_value()) {
                        alpha[2] = az->operator()(xq, yq, zq);
                    }
                }
                //
                // Get the inverse of the Jacobian
                //
                double detJ = 1.0;
                double *__restrict J = &pointJac[dim];
                InverseInPlace<dim>(J, detJ);
                //
                double *__restrict GradN = &NandGradN[nNodes];
                GradPhi.fill(0);
                for (int jn = 0; jn < nNodes; ++jn) {
                    for (int in = 0; in < dim; ++in) {
                        for (int kn = 0; kn < dim; ++kn) {
                            GradPhi[in + jn * dim] += J[in + kn * dim] * GradN[jn + kn * nNodes];
                        }
                    }
                }
                //
                for (int jn = 0; jn < nNodes; ++jn) {
                    for (int in = 0; in <= jn; ++in) {
                        for (int kn = 0; kn < dim; ++kn) {
                            kele[in + jn * nNodes] +=
                                    GradPhi[kn + in * dim] * alpha[kn] * GradPhi[kn + jn * dim] * weight[iq] * detJ;
                        }
                    }
                }
                // Symmetrize the matrix
                for (int jn = 0; jn < nNodes; ++jn) {
                    for (int in = jn + 1; in < nNodes; ++in) {
                        kele[in + jn * nNodes] = kele[jn + in * nNodes];
                    }
                }
                //
                double fq = 0.0;
                if (f.has_value()) {
                    fq = f->operator()(xq, yq, zq);
                }
                for (int in = 0; in < nNodes; ++in) {
                    rele[in] += fq * NandGradN[in] * weight[iq] * detJ;
                }
            }
        }

    };

} // namespace IMSI

//
// Definition of functions
//

namespace IMSI {

    template<typename Device>
    void ScaledLaplacian::GetLinearSystem(Kokkos::View<double *, Device> rhs,
                                          Kokkos::View<size_t *, Device> matRowPtr,
                                          Kokkos::View<int *, Device> matColIdx,
                                          Kokkos::View<double *, Device> matValues) const {

        ///typedef Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space> device_type;

        size_t maxNumDofsPerEle = 0;
        Kokkos::parallel_reduce("MaxDofsPerEle",
                                Kokkos::RangePolicy<Device>(0, meshInfo.mesh.NumberCells()),
                                KOKKOS_LAMBDA(const int &i, size_t &nMax) {
                                    nMax = std::max<size_t>(nMax, size(meshInfo.mesh.NodeList(i)));
                                }, Kokkos::Max<size_t>(maxNumDofsPerEle));

        auto const &c2e = meshInfo.c2e;
        auto const sdim = meshInfo.mesh.GetSpatialDimension();

        for (int ic = 0; ic < c2e.numRows(); ++ic) {
            auto const eleList = c2e.rowConst(ic);
            Kokkos::parallel_for(
                    Kokkos::TeamPolicy<Device>(Kokkos::num_threads(), 1),
                    KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Device>::member_type &team) {
                        auto const ratio = (eleList.length + team.league_size() - 1) / team.league_size();
                        auto const kptr = team.team_rank() + ratio * team.league_rank();
                        std::vector<double> rele(maxNumDofsPerEle);
                        std::vector<double> kele(maxNumDofsPerEle * maxNumDofsPerEle);
                        for (int jj = 0; jj < ratio; ++jj) {
                            auto const ik = kptr + jj;
                            ///
                            //if (ik < eleList.length) {
                            //    printf("COLOR %i Greetings from team %i of league %i --- ratio %i length %i PTR %i \n",
                            //           ic,
                            //           team.team_rank(), team.league_rank(), ratio, eleList.length, ik);
                            //} else {
                            //    printf("COLOR %i Greetings from team %i of league %i --- PTR %i --- SKIPPED \n",
                            //           ic,
                            //           team.team_rank(), team.league_rank(), ik);
                            //}
                            if (ik < eleList.length) {
                                auto const eleID = eleList(ik);
                                auto nodeList = meshInfo.mesh.NodeList(eleID);
                                //
                                // Element type for eleID
                                //
                                switch (meshInfo.mesh.GetCellType(eleID)) {
                                    default:
                                    case ElementType::Q1: {
                                        switch (sdim) {
                                            case 1: {
                                                fe1DQ1 element;
                                                this->ElementaryDataLagrangeFE<1, 2, fe1DQ1>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                            default:
                                            case 2: {
                                                fe2DQ1 element;
                                                this->ElementaryDataLagrangeFE<2, 4, fe2DQ1>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                            case 3: {
                                                fe3DQ1 element;
                                                this->ElementaryDataLagrangeFE<3, 8, fe3DQ1>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                    case ElementType::Q2: {
                                        switch (sdim) {
                                            case 1: {
                                                fe1DQ2 element;
                                                this->ElementaryDataLagrangeFE<1, 3, fe1DQ2>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                            default:
                                            case 2: {
                                                fe2DQ2 element;
                                                this->ElementaryDataLagrangeFE<2, 9, fe2DQ2>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                            case 3: {
                                                fe3DQ2 element;
                                                this->ElementaryDataLagrangeFE<3, 27, fe3DQ2>(element, nodeList, &rele[0],
                                                                                             &kele[0]);
                                                break;
                                            }
                                        }
                                    }
                                }
                                //
                                for (int in = 0; in < size(nodeList); ++in) {
                                    rhs(nodeList[in]) += rele[in];
                                }
                                //
                                for (int in = 0; in < size(nodeList); ++in) {
                                    auto const irow = nodeList[in];
                                    auto const colBegin = &matColIdx(matRowPtr(irow));
                                    auto const colEnd = &matColIdx(matRowPtr(irow + 1));
                                    for (int jn = 0; jn < size(nodeList); ++jn) {
                                        auto const pos = std::lower_bound(colBegin, colEnd, nodeList[jn]) - colBegin;
                                        matValues(matRowPtr(irow) + pos) += kele[in + jn * size(nodeList)];
                                    }
                                }
                            }
                        }
                    }
            );
            Kokkos::fence();
        }
    }

}
