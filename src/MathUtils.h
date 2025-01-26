#pragma once

#include <array>

namespace IMSI {

    /// \brief Inverse N x N matrix in-place and compute determinant
    template<int N, typename Scalar>
    void InverseInPlace(Scalar *__restrict J, Scalar &Jac) {
        if constexpr (N == 1) {
            Jac = J[0];
            J[0] = Scalar(1) / J[0];
        }
        else if constexpr (N == 2) {
            auto const J00 = J[0];
            Jac = J00 * J[3] - J[1] * J[2];
            auto const invJac = Scalar(1) / Jac;
            auto const tmp = J[0];
            //
            J[0] = J[3];
            J[3] = tmp;
            //
            J[1] = -J[1];
            J[2] = -J[2];
            for (int i = 0; i < 4; ++i) {
                J[i] *= invJac;
            }
        }
        else if constexpr (N == 3) {
            auto const iJ33 = J[0]*J[4] - J[1]*J[3];
            auto const iJ32 = J[3]*J[2] - J[5]*J[0];
            auto const iJ31 = J[1]*J[5] - J[2]*J[4];
            //
            auto const iJ23 = J[6]*J[1] - J[7]*J[0];
            auto const iJ22 = J[0]*J[8] - J[2]*J[6];
            auto const iJ21 = J[7]*J[2] - J[1]*J[8];
            //
            auto const iJ13 = J[3]*J[7] - J[4]*J[6];
            auto const iJ12 = J[6]*J[5] - J[8]*J[3];
            auto const iJ11 = J[4]*J[8] - J[5]*J[7];
            //
            Jac = J[0]*iJ11 + J[3]*iJ21 + J[6]*iJ31;
            //
            auto const invJac = Scalar(1)/Jac;
            J[0] = iJ11; J[1] = iJ21; J[2] = iJ31;
            J[3] = iJ12; J[4] = iJ22; J[5] = iJ32;
            J[6] = iJ13; J[7] = iJ23; J[8] = iJ33;
            //
            for (int i = 0; i < 9; ++i) {
                J[i] *= invJac;
            }
        }
        else {
            std::cerr << " Dimension " << N << " x " << N << " not implemented. " << std::endl;
            exit(-1);
        }
    }

}