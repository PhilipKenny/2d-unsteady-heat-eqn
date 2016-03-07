/***************************************************************************************************
 * Name        : 2D_Unsteady_Diffusion.cpp
 * Author      : A. Emre Ongut
 * Version     : 2.0
 * Copyright   : See the copyright notice in the README file.
 * Description : This code is designed for educational purposes.
 *               - For now, everything is for 2D linear triangular elements
 *               - uses Gauss quadrature integration with 7 points
 *               - supports constant Drichlet and Neumann boundary conditions
 *               - not yet in parallel
***************************************************************************************************/

#include "setting.h"
#include "node.h"
#include "element.h"
#include "solver.h"
#include "postProcessor.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
    int numProcs, myRank;
    /***********************************************************************************************
     *  Pre-processing stage: 
     * Read input files, allocate and initialize arrays, etc.
     **********************************************************************************************/
    MPI_Init(&argc, &argv);
    double time = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    setting stngs;              // A setting object is createad.
    stngs.readSettings();       // settings.in file and minf file are read.
    stngs.setProcs(numProcs);
    stngs.setRank(myRank);
     
     
    /**********************************************************************************************/
    size_t nnc, mnc;
    mnc = nnc = (stngs.getNn() - 1)/numProcs + 1;
    if ((myRank + 1)*mnc > stngs.getNn())
        nnc = stngs.getNn() - mnc*myRank; 
    if (mnc*myRank > stngs.getNn())
        nnc = 0;
   
    // printf("Proc # %d from # %d procs with nnLocal = ", myRank, numProcs);
    // cout << nnLocal << "\n";
    node nodes(nnc); // We create a node object
    nodes.createNodes(&stngs); // and read nodal data from file.

    /**********************************************************************************************/
    size_t nec, mec;
    mec = nec = (stngs.getNe() - 1)/numProcs + 1;
    if ((myRank + 1)*mec > stngs.getNe())
        nec = stngs.getNe() - mec*myRank;
    if (mec*myRank > stngs.getNe())
        nec = 0;    

    element elems(nec, nenTri, nefTri, 7);     // We create an element object.
    // printf("Before\n");
    elems.createElements(&stngs, &nodes);                // and read element data from file
    // printf("After\n");
    /**********************************************************************************************/
    // Moreover solution independent calculations are made: Jacobian determinant, shape funct.
    // derivatives as well as element matrices are calculated.

    
    
   
    /***********************************************************************************************
     * Solution stage: 
     * Element matrices are calculated, boundary conditions are applied and system is solved.
     **********************************************************************************************/
    femSolver solver(&stngs, &nodes, &elems);   // We create an femSolver object.
    solver.solve();                             // and solve the equation system

    time = MPI_Wtime() - time;
    MPI_Allreduce(MPI_IN_PLACE, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (stngs.getRank() == 0) cout << time << endl;
    /***********************************************************************************************
     * Post-Processing Stage
     * Create the output files / visualize the reuslts.
     **********************************************************************************************/
//    if(!myRank)
//    {
        postProcessor  postP;
        postP.postProcessorControl(&stngs, &nodes, &elems);
//    }
    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "# " << myRank << "Ciao :)" << endl; 

    MPI_Finalize();


    return 0;
}
