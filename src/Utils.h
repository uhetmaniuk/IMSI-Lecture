#pragma once

namespace IMSI {

    class Mesh;

    void OutputToGMSH
            (
                    const char* fileName,
                    const Mesh& grid,
                    double *p,
                    int numDofs
            );

}
