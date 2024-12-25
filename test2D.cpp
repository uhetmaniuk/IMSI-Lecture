
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>


#include "GlobalVar.h"
#include "Functions_General.h"


#include "SparseSolver.d/linkfc.h"
#include "Utilities.d/GetTime.h"
#include "SparseSolver.d/NgLDLt.hpp"
#include "SparseSolver.d/SymmetricSparse.hpp"


int main(int argc, char **argv)
{

    if (argc==1)
    {
        std::cout << "\n !! Enter the name of an input file !! \n\n";
        exit(-1);
    }

    std::ifstream fin(argv[1]);
    if (!fin) {
        std::cerr << "\n\n !! Input data file " << argv[1] << " opening failed !! \n\n\n";
        return -1;
    }

    ElementType MyElementType = ElementType::Q1;
    int caseNumber = 0;

    int itmp, jtmp, ktmp, ltmp;
    double epsilon = 1.0;

    fin >> itmp >> jtmp >> ktmp >> ltmp;
    switch (itmp)
    {
        case 1:
            MyElementType = ElementType::Q1;
            std::cout << " --> Solve problem with Q1 finite elements \n";
            break;
        case 2:
            MyElementType = ElementType::Q2;
            std::cout << " --> Solve problem with Q2 finite elements \n";
            break;
	    default:
	    	std::cerr << " --> Element Type Not Implemented \n";
	    	return -1;
    }
    fin >> nbElemX >> nbElemY;
    fin >> caseNumber >> epsilon;
    fin.close();

    RealFunction *alpha = 0, *alpha_rhs = 0, *alpha_dir = 0;
    RealD1Function *Solution = 0;
    //--------------
    switch (caseNumber)
    {
        defaut:
        case 0:
            alpha = new Constant(1.0); //--- Laplacian
            alpha_rhs = new Polynomial01();
            Solution = new Parabola01(); //--- Exact solution
            std::cout << " --> Coefficient = Constant \n";
            std::cout << " --> Solution = x (1 - x) y (1 - y) \n";
            break;
        case 1:
            alpha = new CosFunction01;
            alpha_rhs = new CosFunction02;
            Solution = new CosFunction03;
            std::cout << " --> Coefficient and Solution = Varying Pb 1\n";
            break;
        case 2:
            alpha = new CoefficientHou(1.8, epsilon);
            alpha_rhs = new Constant(-1.0);
            std::cout << " --> Epsilon-Coefficient with eps " << epsilon << "\n";
            break;
        case 3:
            alpha = new CoefficientHou_5p1(1.8, 20.0, epsilon);
            alpha_rhs = new Constant(-1.0);
            //alpha_dir = new Constant(0.0);
            alpha_dir = new Parabola00();
            std::cout << " --> Epsilon-Coefficient (Hou & Wu - 5.1) with eps " << epsilon << "\n";
            break;
        case 4:
            alpha = new CoefficientLLL_6p1(epsilon);
            alpha_rhs = new RHS_LLL_6p1();
            std::cout << " --> Epsilon-Coefficient from (LLL - 2013) with eps " << epsilon << "\n";
            break;
        case 5:
            alpha = new CoefficientHou(1.8, epsilon);
            alpha_rhs = new Constant(1.0);
            std::cout << " --> Epsilon-Coefficient (Hou & Wu - 4.4) with eps " << epsilon << "\n";
            std::cout << " --> Test with Neumann Condition \n";
            break;
        case 6:
            alpha = new CoefficientCore(epsilon);
            alpha_rhs = new Constant(1.0);
            std::cout << " --> Core Coefficient with eps " << epsilon << "\n";
            break;
        case 10:
            alpha = new Constant(1.0); //--- Laplacian
            alpha_dir = new SineDir(epsilon);
            std::cout << " --> Laplacian with Oscillatory Dirichlet\n";
            break;
    }
    //--------------

    //
    // Build the mesh (nodes, topology)
    //

    int numNode, numEle, numDofs, locElemSize;

    numEle = nbElemX * nbElemY;

    switch (MyElementType)
    {
        default:
        case ElementType::Q1:
            numNode = (nbElemX + 1) * (nbElemY + 1);
            numDofs = numNode;
            locElemSize = 4;
            break;
        case ElementType::Q2:
            numNode = (2 * nbElemX + 1) * (2 * nbElemY + 1);
            numDofs = numNode;
            locElemSize = 9;
            break;
    }

    std::vector<double> coord_tmp(numNode * 2);
    std::vector<double *> nodeF(numNode);
    for (int ii = 0; ii < numNode; ++ii)
        nodeF[ii] = &coord_tmp[2*ii];

    std::vector<int> elist(numEle * (locElemSize + 2));

    std::vector<int> etmp(locElemSize * numEle);
    int **elemF = new int*[numEle];
    for (int ii = 0; ii < numEle; ++ii)
        elemF[ii] = &etmp[ii * locElemSize];

    std::vector<int> dofPerEle(numEle);

    GetRectangularMesh(0.0, 1.0, 0.0, 1.0,
    		nodeF, elemF,
    		nbElemX, nbElemY, &dofPerEle[0], MyElementType);

    std::cout << " --> End for building rectangular fine mesh\n";

    //
    // Compute the mesh connectivity
    //

    std::vector<int> rowbeg(numDofs + 1);
    std::vector<int> colidx;

    int numNNZ = 0;

    GetSparseMatrixGraph(elemF, numEle, &dofPerEle[0], numDofs,
    		rowbeg, colidx, numNNZ);

    std::cout << " --> End for building sparse graph\n";

    //
    // Build the system matrix
    //

    std::vector<double> val_K(numNNZ, 0.0);
    std::vector<double> val_M;

    double tAssemble = 0.0;
    std::vector<double> tMake(16, 0.0);

    //--- Assemble the stiffness matrix
    tAssemble -= getTime();
    GetMatrixData(nodeF, elemF, numEle,
                  rowbeg, colidx, numNNZ, val_K, val_M,
                  MyElementType, alpha, &tMake[0]);
    tAssemble += getTime();

    //--- Matrix KK stores the stiffness matrix for the discretization
    //--- without any boundary condition
    //--- Q1, Q2: KK = K_{fem}
    //--- MFEM_L: KK = Phi^T K_{fine} Phi
    SymmetricSparse<double> KK(numDofs, numNNZ, &rowbeg[0], &colidx[0], &val_K[0]);

    std::cout << " --> End for building KK\n";

    //
    // Build the projection of Dirichlet condition
    //

    std::vector<double> dirF(numDofs, 0.0);
    std::vector<int> globalToFree(numDofs, -1);
    int numFree = 0;

    double timeDir = -getTime();
    if (caseNumber == 5) {

        if (alpha_dir)
        	exit(-1);

        double htol = 1.0e-03 / (nbElemX + nbElemY);
        for (int ii = 0; ii < numNode; ++ii) {
            double xx = nodeF[ii][0], yy = nodeF[ii][1];
            if ((std::abs(xx) < htol) || (std::abs(xx - 1.0) < htol)
                || (std::abs(yy - 1.0) < htol))
                continue;
            globalToFree[ii] = numFree++;
        } // for (int ii = 0; ii < numNode; ++ii)

    }
    else {
        GetDirichletBC(nodeF, elemF, MyElementType, globalToFree,
                       numFree, dirF, alpha_dir);
    }
    timeDir += getTime();

    std::cout << " --> End for setting the Dirichlet condition right hand side\n";

    //
    // Build the right hand side
    //

    std::vector<double> rhsF(numDofs, 0.0);

    if (caseNumber == 5) {
        if (MyElementType == ElementType::Q1) {
            const double hh = double(1.0)/double(nbElemX);
            const double xi = -1.0/sqrt(3.0);
            for (int ii = 0; ii < nbElemX; ++ii) {
                double myx = (ii + 0.5 + 0.5*xi) * hh;
                rhsF[ii] += 0.5 * hh * sin(M_PI*myx) * (-xi + 1)/2;
                rhsF[ii+1] += 0.5 * hh * sin(M_PI*myx) * (xi + 1)/2;
                myx = (ii + 0.5 - 0.5*xi) * hh;
                rhsF[ii] += 0.5 * hh * sin(M_PI*myx) * (xi + 1)/2;
                rhsF[ii+1] += 0.5 * hh * sin(M_PI*myx) * (-xi + 1)/2;
            }
        }
        else if (MyElementType == ElementType::Q2) {
            //
            std::cout << "\n !!! Q2 discretization not implemented !!! \n\n";
            assert(0 > 1);
            //
            const double hh = 1.0/double(nbElemX);
            for (int ii = 0; ii < nbElemX; ++ii) {
                //rhsF[ii] += 0.5*hh;
                //rhsF[ii+1] += 0.5*hh;
                //rhsF[ii+1] += hh/6.0;
            }
        }
    }
    else {
        GetRightHandSideData(nodeF, elemF, numEle, rhsF, MyElementType, alpha_rhs);
    }

    std::cout << " --> End for building the right hand side\n";

    std::vector<int> freeToGlobal(numFree);
    for (int id = 0; id < numDofs; ++id)
    {
        if (globalToFree[id] > -1)
            freeToGlobal[globalToFree[id]] = id;
    }

    //
    // Update the right hand side with a Dirichlet condition
    // K * [0; g_b] = [K_{ib} * g_{b}; K_{bb} * g_{b}]
    //

    if (alpha_dir)
    {
        std::vector<double> Kd(numDofs, 0.0);
        KK.Apply(1, &dirF[0], &Kd[0]);
        for (int id = 0; id < numDofs; ++id) {
            if (globalToFree[id] > -1)
                rhsF[id] -= Kd[id];
        }
    }

    //
    // Solve the resulting algebraic system
    //

    std::vector<double> u(numDofs, 0.0);

    double tFactor = 0.0, tSolve = 0.0;
    {
        SymmetricSparse<double> *Ktmp = KK.SymmetricExtract(numFree, &freeToGlobal[0]);
        std::cout << " --> End for building the system matrix\n";
        tFactor -= getTime();
        Ktmp->factor();
        tFactor += getTime();
        std::cout << " --> End for factoring K\n";
        //
        NgLDLt<double, int> mysolver(*Ktmp);
        mysolver.factorize(*Ktmp);
        //----------
        std::vector<double> rhs_u(numFree, 0.0);
        for (int ii = 0; ii < numDofs; ++ii)
        {
            if (globalToFree[ii] > -1)
                rhs_u[globalToFree[ii]] = rhsF[ii];
        }
        //----------
        std::vector<double> u_u(numFree, 0.0);
        tSolve -= getTime();
        Ktmp->Solve(1, &rhs_u[0], &u_u[0]);
        tSolve += getTime();
        //----------
        std::vector<double> vv(numFree, 0.0);
        rhs_u.assign(numFree, 0.0);
        for (int ii = 0; ii < numDofs; ++ii)
        {
            if (globalToFree[ii] > -1)
                rhs_u[globalToFree[ii]] = rhsF[ii];
        }
        mysolver.Apply(1, &rhs_u[0], &vv[0]);
        for (int ii = 0; ii < 13; ++ii) {
            std::cout << ii << " " << vv[ii] << " " << u_u[ii] << "\n";
        }
        for (int ii = numFree-19; ii < numFree; ++ii) {
            std::cout << ii << " " << vv[ii] << " " << u_u[ii] << "\n";
        }
        //----------
        for (int ii = 0; ii < numDofs; ++ii)
        {
            if (globalToFree[ii] > -1)
                u[ii] = u_u[globalToFree[ii]];
            else {
                if (alpha_dir) {
                    u[ii] = dirF[ii];
                }
            }
        } // for (int ii = 0; ii < numDofs; ++ii)
        //----
        delete Ktmp;
    }

    std::cout << " --> End for solving\n";

    //
    // Print information on the screen
    //

    std::cout << "------------------------------------------------\n";
    std::cout << "Nb of x-elements             = " << nbElemX << std::endl;
    std::cout << "Nb of y-elements             = " << nbElemY << std::endl;
    std::cout << "Nb of nodes                  = " << numNode << std::endl;
    std::cout << "Nb of degrees of freedom     = " << numDofs << std::endl;
    std::cout << "X-meshsize                   = " << 1.0/nbElemX << std::endl;
    std::cout << "Y-meshsize                   = " << 1.0/nbElemY << std::endl;
    std::cout << "------------------------------------------------\n";
    std::cout << "Assembly (s)                 = " << tAssemble/1000.0 << std::endl;
    std::cout << "Dirichlet Right Hand Side (s)= " << timeDir/1000.0 << std::endl;
    std::cout << "Factorization (s)            = " << tFactor/1000.0 << std::endl;
    std::cout << "Solve (s)                    = " << tSolve/1000.0 << std::endl;
    std::cout << "------------------------------------------------\n";

    //
    // Compute the nodal values of the exact solution
    //

    if (Solution)
    {

        double domainL2Norm = 0.0, domainL2Err = 0.0;
        double domainH1Norm = 0.0, domainH1Err = 0.0;
        double difMax = 0.0;
        double refMax = 0.0;
        int nodeMax = 0;

        //
        //--- Compute the L^max norm of error on the nodes
        //
        for (int j=0; j<numNode; ++j)
        {
            std::vector<double> point(2);
            point[0] = nodeF[j][0]; point[1] = nodeF[j][1];
            double uExact = Solution->Apply(point);
            //----
            double tmp = std::abs(uExact);
            refMax = (tmp > refMax) ? tmp : refMax;
            //---
            tmp = std::abs(uExact - u[j]);
            if (tmp > difMax)
            {
                difMax = tmp;
                nodeMax = j;
            }
        } // for (int j=0; j<numNode; ++j)

        //
        // Compute the L^2 and H^1 error over the domain
        //

        GetErrorNorms(nodeF, elemF, numEle, &dofPerEle[0], MyElementType,
                      u, Solution,
                      domainL2Err, domainL2Norm, domainH1Err, domainH1Norm);

        std::cout.precision(4);
        std::cout << "Error absolute L^2 in domain = " << domainL2Err << std::endl;
        std::cout << "Error relative L^2 in domain = " << domainL2Err/domainL2Norm << std::endl;
        std::cout << "------------------------------------------------\n";
        std::cout << "Error absolute H^1 in domain = " << domainH1Err << std::endl;
        std::cout << "Error relative H^1 in domain = " << domainH1Err/domainH1Norm << std::endl;
        std::cout << "------------------------------------------------\n";
        std::cout << "Max Error absolute on nodes  = " << difMax << std::endl;
        std::cout << "Max Error relative on nodes  = " << difMax/refMax << std::endl;
        std::cout << "Node with MAX error          = (" << nodeF[nodeMax][0];
        std::cout << ", " << nodeF[nodeMax][1] << ")\n";
        std::cout << "------------------------------------------------\n";

    } // if (Solution)

    //
    //--- Output solution to GMSH
    //

    double energy = 0.0;

        //
        //--- Treat the case of Q1, Q2 discretization
        //

        OutputMesh("Test2D.msh", numNode, nodeF, numEle, elemF, MyElementType, &u[0]);

        if (alpha_dir) {
            std::vector<double> utilde(numDofs), Kutilde(numDofs);
            for (int ii = 0; ii < numDofs; ++ii) {
                std::vector<double> mynode(2);
                mynode[0] = nodeF[ii][0]; mynode[1] = nodeF[ii][1];
                double g_ii = alpha_dir->Apply(mynode);
                utilde[ii] = u[ii] - g_ii;
            }
            //
            KK.Apply(1, &utilde[0], &Kutilde[0]);
            //
            for (int ii = 0; ii < numDofs; ++ii)
                rhsF[ii] = 0.0;
            GetRightHandSideData(nodeF, elemF, numEle, rhsF, MyElementType, alpha_rhs);
            //
            for (int ii = 0; ii < numDofs; ++ii) {
                std::vector<double> mynode(2);
                mynode[0] = nodeF[ii][0]; mynode[1] = nodeF[ii][1];
                double g_ii = alpha_dir->Apply(mynode);
                energy += utilde[ii] * (0.5 * Kutilde[ii] - rhsF[ii]) + g_ii * Kutilde[ii];
            } // for (int ii = 0; ii < numDofs; ++ii)
        }
        else {
            std::vector<double> Ku(numDofs);
            KK.Apply(1, &u[0], &Ku[0]);
            for (int ii = 0; ii < numDofs; ++ii)
                energy += u[ii] * (0.5 * Ku[ii] - rhsF[ii]);
        }

    std::cout.precision(16);
    std::cout << "Energy                       = " << energy << std::endl;
    std::cout << "------------------------------------------------\n";

    //--- Delete arrays

    delete[] elemF;

    return 0;

}

