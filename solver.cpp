#include "solver.h"
#include "mpi.h"
#include <omp.h>

/***************************************************************************************************
 * void femSolver::solve()
 * Driver routine for the femSolver
***************************************************************************************************/
void femSolver::solve()
{ 
    postProcessor  postP;
    nodes->exchangeTemp();
    unsigned t = 0;
    unsigned max_t = stngs->getNIter();
    double tol = stngs->getTolerance();
    while(t++ < max_t && !explicitSolver(tol))
//        if(t%2500 == 0)postP.postProcessorControl(stngs, nodes, elems)
;


//    double time = (double)clock()/CLOCKS_PER_SEC;
//    MPI_Allreduce(MPI_IN_PLACE, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//    if (stngs->getRank() == 0) cout << time << endl;
//    cout << endl << "Total execution time is " << fixed << (float)clock()/CLOCKS_PER_SEC 
    //             << " seconds." << endl;

    return;
}


/***************************************************************************************************
 * void femSolver::explicitSolver()
 * This routine solves the equation system using explicit integration.
***************************************************************************************************/
bool femSolver::explicitSolver(const double& tolerance)
{
    size_t totNodes = nodes->getNn() + nodes->getAdditionalNodes();
//if(nodes->getRank()==1)for(int i = 0; i< totNodes; ++i) cout<<i<<"\t"<<nodes->getT(i) <<endl;
//cout<<endl;
    elems -> calculateRHS(nodes, stngs);
    double *temp = new double[totNodes]();

//if(nodes->getRank()==1)for(int i = 0; i< totNodes; ++i) cout<<nodes->getT(i) <<endl;

//if(nodes->getRank()==1)for(int i = 0; i< totNodes; ++i) cout<<temp[i] <<endl;
#pragma omp parallel for
    for(size_t ie = 0; ie < elems->getNe(); ++ie)
    {
        for(int ien = 0; ien < elems -> getNen(); ++ien)
        {
//if(nodes->getRank()==1)cout<<elems->getConn(ie,ien)<<endl;
            double* RHS = elems->getRHSptr(ie);
            size_t localNode = ie*(elems -> getNen()) + ien;
            size_t globalNode = elems -> getConn(ie, ien);
	    temp[globalNode] += RHS[ien] / nodes -> getMass(globalNode);
//            if(nodes->getRank()==1)cout << globalNode << "\t" << temp[globalNode] <<endl;
	}
   } 
#pragma omp parallel for
    for (size_t in = 0; in < totNodes; ++in) nodes->setT(in, 0);

    int isConverged = true;
//    for (size_t in =0; in < totNodes; ++in)
//    {
#pragma omp parallel for
    for(size_t ie = 0; ie < elems->getNe(); ++ie)
    {
        for(int ien = 0; ien < elems -> getNen(); ++ien)
        {
        size_t globalNode = elems->getConn(ie,ien);
        if(fabs(temp[globalNode] - nodes->getT(globalNode)) > tolerance) isConverged = false;
	nodes -> setT(globalNode, temp[globalNode]); 
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &isConverged, 1, MPI_INT, MPI_BAND, MPI_COMM_WORLD);
    nodes->accumulateTemp();
    nodes->exchangeTemp();
//    for(int i = 0; i < totNodes; ++i) cout<<nodes->getRank()<<"\t"<<i<<"\t"<<temp[i] <<endl;
    delete[] temp;
    return isConverged;
}


