#include <iostream>
#include <vector>

#include "QuadratureRule.h"

namespace IMSI {

    void getQuadrature
            (RuleType type,
             int sdim,
             int order,
             int &length,
             std::vector<double> &weights,
             std::vector<double> &xi,
             std::vector<double> &eta,
             std::vector<double> &zeta
            ) {

        constexpr int padding = 4;

        if ((sdim == 1) && (type == RuleType::Gauss)) {

            length = order;
            weights.resize(((length + padding - 1) / padding) * padding);
            xi.resize(weights.size());
            eta.assign(weights.size(), 0);
            zeta.assign(weights.size(), 0);

            switch (order) {
                default:
                case 1: {
                    xi[0] = 0.0;
                    weights[0] = 2.0;
                    break;
                }
                case 2: {
                    xi[0] = -1.0 / sqrt(3.0);
                    xi[1] = -xi[0];
                    //----
                    weights[0] = 1.0;
                    weights[1] = 1.0;
                    break;
                }
                case 3: {
                    xi[0] = -sqrt(0.6);
                    xi[1] = 0.0;
                    xi[2] = -xi[0];
                    weights[0] = 5.0 / 9.0;
                    weights[1] = 8.0 / 9.0;
                    weights[2] = weights[0];
                    break;
                }
                case 4: {
                    xi[0] = -0.861136311594053;
                    xi[1] = -0.339981043584856;
                    xi[2] = -xi[1];
                    xi[3] = -xi[0];
                    weights[0] = 0.347854845137454;
                    weights[1] = 0.652145154862546;
                    weights[2] = weights[1];
                    weights[3] = weights[0];
                    break;
                }
                case 5: {
                    xi[0] = -0.906179845938664;
                    xi[1] = -0.538469310105683;
                    xi[2] = 0.0;
                    xi[3] = -xi[1];
                    xi[4] = -xi[0];
                    weights[0] = 0.236926885056189;
                    weights[1] = 0.478628670499366;
                    weights[2] = 0.568888888888889;
                    weights[3] = weights[1];
                    weights[4] = weights[0];
                    break;
                }
                case 6: {
                    xi[0] = -0.932469514203152;
                    xi[1] = -0.661209386466265;
                    xi[2] = -0.238619186083197;
                    xi[3] = -xi[2];
                    xi[4] = -xi[1];
                    xi[5] = -xi[0];
                    weights[0] = 0.171324492379170;
                    weights[1] = 0.360761573048139;
                    weights[2] = 0.467913934572691;
                    weights[3] = weights[2];
                    weights[4] = weights[1];
                    weights[5] = weights[0];
                    break;
                }
                case 7: {
                    xi[0] = -0.949107912342759;
                    xi[1] = -0.741531185599394;
                    xi[2] = -0.405845151377397;
                    xi[3] = 0.0;
                    xi[4] = -xi[2];
                    xi[5] = -xi[1];
                    xi[6] = -xi[0];
                    weights[0] = 0.129484966168870;
                    weights[1] = 0.279705391489277;
                    weights[2] = 0.381830050505119;
                    weights[3] = 0.417959183673469;
                    weights[4] = weights[2];
                    weights[5] = weights[1];
                    weights[6] = weights[0];
                    break;
                }
                case 8: {
                    xi[0] = -0.960289856497536;
                    xi[1] = -0.796666477413627;
                    xi[2] = -0.525532409916329;
                    xi[3] = -0.183434642495650;
                    xi[4] = -xi[3];
                    xi[5] = -xi[2];
                    xi[6] = -xi[1];
                    xi[7] = -xi[0];
                    weights[0] = 0.101228536290376;
                    weights[1] = 0.222381034453374;
                    weights[2] = 0.313706645877887;
                    weights[3] = 0.362683783378362;
                    weights[4] = weights[3];
                    weights[5] = weights[2];
                    weights[6] = weights[1];
                    weights[7] = weights[0];
                    break;
                }
                case 9: {
                    xi[0] = -0.968160239507626;
                    xi[1] = -0.836031107326636;
                    xi[2] = -0.613371432700590;
                    xi[3] = -0.324253423403809;
                    xi[4] = 0.0;
                    xi[5] = -xi[3];
                    xi[6] = -xi[2];
                    xi[7] = -xi[1];
                    xi[8] = -xi[0];
                    weights[0] = 0.081274388361574;
                    weights[1] = 0.180648160694857;
                    weights[2] = 0.260610696402935;
                    weights[3] = 0.312347077040003;
                    weights[4] = 0.330239355001260;
                    weights[5] = weights[3];
                    weights[6] = weights[2];
                    weights[7] = weights[1];
                    weights[8] = weights[0];
                    break;
                }
                case 10: {
                    xi[0] = -0.973906528517172;
                    xi[1] = -0.865063366688985;
                    xi[2] = -0.679409568299024;
                    xi[3] = -0.433395394129247;
                    xi[4] = -0.148874338981631;
                    xi[5] = -xi[4];
                    xi[6] = -xi[3];
                    xi[7] = -xi[2];
                    xi[8] = -xi[1];
                    xi[9] = -xi[0];
                    weights[0] = 0.066671344308688;
                    weights[1] = 0.149451349150581;
                    weights[2] = 0.219086362515982;
                    weights[3] = 0.269266719309996;
                    weights[4] = 0.295524224714753;
                    weights[5] = weights[4];
                    weights[6] = weights[3];
                    weights[7] = weights[2];
                    weights[8] = weights[1];
                    weights[9] = weights[0];
                    break;
                }
            } // switch(order)

        } else if ((sdim == 1) && (type == RuleType::GaussLobatto)) {

            int myOrder = (order > 1) ? order : 2;
            length = myOrder;
            weights.resize(((length + padding - 1) / padding) * padding);
            xi.resize(weights.size());
            eta.assign(weights.size(), 0);
            zeta.assign(weights.size(), 0);

            switch (myOrder) {
                default:
                case 2: {
                    xi[0] = -1.0;
                    xi[1] = 1.0;
                    weights[0] = 1.0;
                    weights[1] = 1.0;
                    break;
                }
                case 3: {
                    xi[0] = -1.0;
                    xi[1] = 0.0;
                    xi[2] = 1.0;
                    weights[0] = 1.0 / 3.0;
                    weights[1] = 4.0 * weights[0];
                    weights[2] = weights[0];
                    break;
                }
            }

        }

        if (sdim == 1)
            return;

        if ((type == RuleType::Gauss) || (type == RuleType::GaussLobatto)) {

            std::vector<double> p1d, w1d;
            int length1D = 0;
            getQuadrature(type, 1, order, length1D, w1d, p1d, eta, zeta);

            if (sdim == 2) {
                length = length1D * length1D;
                size_t const myLen = ((length + padding - 1) / padding) * padding;
                weights.assign(myLen, 0);
                xi.assign(weights.size(), 0);
                eta.assign(weights.size(), 0);
                zeta.assign(weights.size(), 0);
                size_t counter = 0;
                for (size_t jq = 0; jq < length1D; ++jq) {
                    for (size_t iq = 0; iq < length1D; ++iq) {
                        weights[counter] = w1d[iq] * w1d[jq];
                        xi[counter] = p1d[iq];
                        eta[counter] = p1d[jq];
                        counter += 1;
                    }
                }
            } else if (sdim == 3) {
                length = w1d.size() * w1d.size() * w1d.size();
                size_t myLen = ((length + padding - 1) / padding) * padding;
                weights.assign(myLen, 0);
                xi.assign(weights.size(), 0);
                eta.assign(weights.size(), 0);
                zeta.assign(weights.size(), 0);
                size_t counter = 0;
                for (size_t kq = 0; kq < length1D; ++kq) {
                    for (size_t jq = 0; jq < length1D; ++jq) {
                        for (size_t iq = 0; iq < length1D; ++iq) {
                            weights[counter] = w1d[iq] * w1d[jq] * w1d[kq];
                            xi[counter] = p1d[iq];
                            eta[counter] = p1d[jq];
                            zeta[counter] = p1d[kq];
                            counter += 1;
                        }
                    }
                }
            } else {
                std::cout << "\n !!! getQuadrature >> ";
                std::cout << " Dimension " << sdim
                          << " is not implemented !!!\n\n";
                assert((sdim == 2) || (sdim == 3));
            }
            return;

        }

    }

}