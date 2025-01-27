#include "Element.h"
#include "fe1DQ1.h"
#include "fe1DQ2.h"
#include "fe2DQ1.h"
#include "fe2DQ2.h"
#include "fe3DQ1.h"
#include "fe3DQ2.h"
#include "Mesh.h"
#include "QuadratureRule.h"
#include "Utils.h"

#include "Kokkos_Core.hpp"

#include <array>
#include <fstream>
#include <ostream>

namespace IMSI {

    void MapDegreesOfFreedom(const std::vector<int>& bdyNodes,
                             std::vector<int>& globalToFree,
                             std::vector<int>& freeToGlobal) {
        globalToFree.assign(size(globalToFree), 0);
        for (auto ib : bdyNodes) {
            globalToFree[ib] = -1;
        }
        int icount = 0;
        for (int i = 0; i < size(globalToFree); ++i) {
            if (globalToFree[i] == -1) {
                continue;
            }
            globalToFree[i] = icount;
            freeToGlobal[icount] = i;
            icount += 1;
        }
    }

    /*
    namespace {

        template<int dim, int nNodes, typename ElementClass>
        void ElementaryNorms(ElementClass &element,
                                      const std::vector<int> &nodeList,
                             double *u,
                             const std::optional< std::function<double(double, double, double)> >& solution,
                             double values[4]) {

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
                double solution_q = (solution.has_value()) ? solution->operator()(xq, yq, zq) : 0.0;
                double uq = 0.0;
                for (int in = 0; in < nNodes; ++in) {
                    uq += u[nodeList[in]] * NandGradN[in];
                }
                //
                values[0] += (solution_q - uq) * (solution_q - uq) * weight[iq] * detJ;
                values[1] += solution_q * solution_q * weight[iq] * detJ;
                values[2] += 0 * weight[iq] * detJ;
                values[3] += 0 * weight[iq] * detJ;
            }
        }
    }

    std::tuple<double, double, double, double> GetErrorNorms( const Mesh& grid,
                                                              double *u,
                                                              const std::optional< std::function<double(double, double, double)> >& solution) {

        struct NormReduct {
            double values[4];
            void operator+=(NormReduct const& other) {
                for (int i = 0; i < 4; ++i) {
                    values[i] += other.values[i];
                }
            }
            // In multi-threaded environments, variables shared between threads can be modified by other threads at any time.
            // Using volatile ensures that the compiler doesn't make incorrect assumptions about the variable's
            // value based on its own thread's execution.
            void operator+=(NormReduct const volatile& other) volatile {
                for (int i = 0; i < 4; ++i) {
                    values[i] += other.values[i];
                }
            }
        };
        NormReduct estimate = {{0, 0, 0, 0}};

        Kokkos::parallel_reduce(grid.NumberCells(), KOKKOS_LAMBDA(int iE, NormReduct &valtmp) {

            auto const& nodeList = grid.NodeList(iE);
            //
            // Element type for eleID
            //
            switch (grid.GetCellType(iE)) {
                default:
                case ElementType::Q1: {
                    switch (grid.GetSpatialDimension()) {
                        case 1: {
                            fe1DQ1 element;
                            ElementaryNorms<1, 2, fe1DQ1>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                        default:
                        case 2: {
                            fe2DQ1 element;
                            ElementaryNorms<2, 4, fe2DQ1>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                        case 3: {
                            fe3DQ1 element;
                            ElementaryNorms<3, 8, fe3DQ1>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                    }
                    break;
                }
                case ElementType::Q2: {
                    switch (grid.GetSpatialDimension()) {
                        case 1: {
                            fe1DQ2 element;
                            ElementaryNorms<1, 3, fe1DQ2>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                        default:
                        case 2: {
                            fe2DQ2 element;
                            ElementaryNorms<2, 9, fe2DQ2>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                        case 3: {
                            fe3DQ2 element;
                            ElementaryNorms<3, 27, fe3DQ2>(element, nodeList, u, solution, valtmp.values);
                            break;
                        }
                    }
                }
            }
        }, estimate);

        for (double & value : estimate.values) {
            value = std::sqrt(value);
        }

        return {estimate.values[0], estimate.values[1], estimate.values[2], estimate.values[3]};

    }
*/

    void OutputToGMSH
            (
                    const char* fileName,
                    const Mesh& grid,
                    double *p,
                    int numDofs
            )
    {

        std::ofstream fout(fileName);

        fout << "$MeshFormat" << std::endl;
        fout << "2.0 0 8\n";
        fout << "$EndMeshFormat\n";

        fout << "\n";

        fout << "$Nodes\n";
        fout << grid.NumberVertices() << std::endl;
        for (int iN = 0; iN < grid.NumberVertices(); ++iN) {
            fout << iN+1 << " ";
            auto const coord = grid.GetVertex(iN);
            fout << coord[0] << " " << coord[1] << " " << coord[2];
            fout << std::endl;
        }
        fout << "$EndNodes" << std::endl;

        fout << "$Elements\n";
        fout << grid.NumberCells() << std::endl;
        for (int iE = 0; iE < grid.NumberCells(); ++iE) {
            fout << iE + 1;
            switch (grid.GetCellType(iE)) {
                default:
                case ElementType::Q1: {
                    fout << " 3 0 ";
                    break;
                }
                case ElementType::Q2: {
                    fout << " 10 0 ";
                    break;
                }
            }
            auto const nodeList = grid.NodeList(iE);
            for (auto iN : nodeList) {
                fout << iN + 1 << " ";
            }
            fout << std::endl;
        }
        fout << "$EndElements" << std::endl;

        fout << "\n";

        if (p) {
            fout << "$NodeData\n";
            fout << "1\n";
            fout << "\"Scalar Field\"\n";
            fout << "1\n";
            fout << "0.0\n";
            fout << "3\n";
            fout << "0\n";
            fout << "1\n";

            fout << numDofs << "\n";
            for (int in = 0; in < numDofs; ++in)
            {
                fout << in + 1 << " " << p[in] << std::endl;
            }
            fout << "$EndNodeData\n";
        }

        fout.close();

    }

}
