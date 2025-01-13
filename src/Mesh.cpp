#include "Mesh.h"

#include <ostream>

namespace IMSI {

    std::ostream &operator<<
            (
                    std::ostream &os, const IMSI::Mesh &ref
            ) {

        os << "//\n// NODES\n//\n";

        for (int iN = 0; iN < ref.vertexList.size(); ++iN) {
            os.precision(8);
            os.setf(std::ios::scientific, std::ios::floatfield);
            auto const &coord = ref.vertexList[iN];
            for (int jD = 0; jD < 3; ++jD)
                os << coord[jD] << " ";
            os << std::endl;
        } // for (int iN = 0; iN < ref.d_numVertices; ++iN)

        os << "//\n// MESH CONNECTIVITY\n//\n";

        for (int iE = 0; iE < ref.cellToNode.size(); ++iE) {
            auto const &vertexList = ref.cellToNode[iE];
            for (int jN = 0; jN < vertexList.size(); ++jN)
                os << vertexList[jN] << " ";
            os << std::endl;
        } // for (int iE = 0; iE < ref.d_numCells; ++iE)

        return os;

    }

}


/*

//
// Questions:
// * Do we need a global side list?
// * How to manage the orientation of sides?????
// * Do we need to store sideList in each element?
// * How to create the map sideToCell?
// * How to identify edges in 3D?
// * Do we need to store elemToElem?
//

*/

