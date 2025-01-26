#include "Mesh.h"

#include <ostream>

namespace IMSI {

    Mesh::Mesh(
            int dim,
            std::vector< std::array<double, 3> > const& vertex,
            std::vector<ElementType>&& cList,
            std::vector< std::vector<int> >&& cToN,
            std::vector<int>&& nodeBdry
    ) : sdim(dim), vertex_x(vertex.size()), vertex_y(vertex.size()), vertex_z(vertex.size()),
        cellType(std::move(cList)), cellToNode(std::move(cToN)),
        boundaryNode(std::move(nodeBdry)) {
        size_t count = 0;
        for (auto const& coord : vertex) {
            vertex_x[count] = coord[0];
            vertex_y[count] = coord[1];
            vertex_z[count] = coord[2];
            count += 1;
        }
    }

    std::ostream &operator<<
            (
                    std::ostream &os, const IMSI::Mesh &ref
            ) {

        os << "//\n// NODES\n//\n";

        for (int iN = 0; iN < ref.vertex_x.size(); ++iN) {
            os.precision(8);
            os.setf(std::ios::scientific, std::ios::floatfield);
            os << ref.vertex_x[iN] << " " << ref.vertex_y[iN] << " " << ref.vertex_z[iN];
            os << std::endl;
        } // for (int iN = 0; iN < ref.d_numVertices; ++iN)

        os << "//\n// MESH CONNECTIVITY\n//\n";

        for (const auto & vertexList : ref.cellToNode) {
            for (int jN : vertexList) {
                os << jN << " ";
            }
            os << std::endl;
        }

        return os;

    }

}
