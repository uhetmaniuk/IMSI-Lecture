#include "Mesh.h"
#include "Utils.h"

#include <array>
#include <fstream>
#include <ostream>

namespace IMSI {

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
