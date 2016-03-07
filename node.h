#ifndef NODE_H_
#define NODE_H_

#include "setting.h"
#include "mpi.h"
#include <vector>
#include <algorithm> 

const unsigned numWins = 5;
/***************************************************************************************************
 * GENERAL NODE LEVEL DATA STRUCTURE
 ***************************************************************************************************
 * Here we store the nodal data in a structure of arrays form.
 **************************************************************************************************/
class node
{
    private:
        // PRIVATE VARIABLES
        size_t  nn;     // Number of nodes
        size_t  additionalNodes;
        std::vector<int>    BCtype; // Boundary condition type
        std::vector<double> x;      // x-coordinate array pointer
        std::vector<double> y;      // y-ccordinate array pointer
        std::vector<double> T;      // Temperature array pointer
        std::vector<double> mass;   // nodal mass

        // MPI VARIABLES
        MPI_Win win[numWins];                  // Communication windows
        std::vector< vector<int> > offsets;    // Positions of data on remote process:w

//        std::vector< vector<int> > localOffsets;    // Positions of additional data on local process
/*        std::vector< vector<int> > remoteBC;   // Container for remote node BCs
        std::vector< vector<double> > remoteX; // Container for remote node x-coordinates
        std::vector< vector<double> > remoteY; // Container for remote node y-coordinates
        std::vector< vector<double> > remoteT; // Container for remote node temperature
        std::vector< vector<double> > remoteM; // Container for remote node mass
*/
        MPI_Datatype *remoteDatatypeInt;
        MPI_Datatype *remoteDatatype;
        MPI_Datatype *localDatatypeInt;
        MPI_Datatype *localDatatype;
        size_t maxNodes;
        int myRank;
        int numProcs;
            
        // PRIVATE METHODS
        void readMxyz(setting* settings);
        void readData(setting* settings);  
        void readMxyzPar(setting* settings);
        void readDataPar(setting* settings);
        void assignRemoteStorage();
        void createDatatypes();
        void deleteDatatypes(); 
        void createWindows();
        void freeWindows(int numWins);
        void exchangeData();
        void sortOffsets(); 
 
    protected:

    public:
        // CONSTRUCTOR
        node(size_t argNn)
        {
            nn     = argNn;
            additionalNodes = 0;
            BCtype.resize(nn);
            x.resize(nn);   //   = new double[nn];
            y.resize(nn);   //   = new double[nn];
            T.resize(nn);   //   = new double[nn];
            mass.resize(nn);//   = new double[nn];
        };

        // DESTRUCTOR
        ~node()
        {

           // delete[] BCtype;
           // delete[] x;
           // delete[] y;
           // delete[] T;
           // delete[] mass;
        };

        // GETTERS
        inline double getX(size_t i)      {return x[i];};
        inline double getY(size_t i)      {return y[i];};
        inline double getT(size_t i)      {return T[i];};
        inline double getMass(size_t i)   {return mass[i];};
        inline int    getBCtype(size_t i) {return BCtype[i];};
        inline size_t getNn()             {return nn;};
        inline size_t getAdditionalNodes(){return additionalNodes;};
        inline int    getRank()           {return myRank;};
        inline size_t getMaxNodes()       {return maxNodes;};

	void initMPI();
//        void sortOffsets(); 
        int getRemoteNode(size_t globalNode){return globalNode % maxNodes;};
        int getRemoteRank(size_t globalNode){return globalNode / maxNodes;};
	size_t searchArr(size_t in); // Returns index should not be here
        size_t getExtra(int rank){return offsets[numProcs][rank];};
      /*  
	double getX(size_t in)
        {
          int remoteRank = getRemoteRank(in);
	  int remoteNode = getRemoteNode(in);
          if (remoteRank == myRank)
            return x[remoteNode];
          else
            return remoteX[remoteRank][searchArr(remoteNode)];
        };

	double getY(size_t in)
        {
          if (getRemoteRank(in) == myRank)
            return y[in];
          else
            return remoteY[getRemoteRank(in)][searchArr(in)];
        };

        double getT(size_t in)
        {
          if (getRemoteRank(in) == myRank)
            return T[in];
          else
            return remoteT[getRemoteRank(in)][searchArr(in)];
        };

        double getMass(size_t in)
        {
          if (getRemoteRank(in) == myRank)
            return mass[in];
          else
            return remoteM[getRemoteRank(in)][searchArr(in)];
        };

        int getBCtype(size_t in)
        {
          if (getRemoteRank(in) == myRank)
            return BCtype[in];
          else
            return remoteBC[getRemoteRank(in)][searchArr(in)];
        };
*/

        // SETTERS 
         inline void setMass   (size_t i, double value) {mass[i]    = value;};
         inline void addMass   (size_t i, double value) {mass[i]   += value;};
         inline void setT      (size_t i, double value) {T[i]       = value;};
         inline void setBCtype (size_t i, int    value) {BCtype[i]  = value;};
/*
        void setMass(size_t in, double value)
        {
          if (getRemoteRank(in) == myRank)
            mass[in] = value;
          else
            remoteM[getRemoteRank(in)][searchArr(in)] = value;
        };

        void addMass(size_t in, double value)
        {
          if (getRemoteRank(in) == myRank)
            mass[in] += value;
          else
            remoteM[getRemoteRank(in)][searchArr(in)] += value;
        };

        void setT(size_t in, double value)
        {
          if (getRemoteRank(in) == myRank)
            T[in] = value;
          else
            remoteT[getRemoteRank(in)][searchArr(in)] = value;
        };

        void setBCtype(size_t in, double value)
        {
          if (getRemoteRank(in) == myRank)
            BCtype[in] = value;
          else
            remoteBC[getRemoteRank(in)][searchArr(in)] = value;
        };

*/
	// @ Klim
	inline double  getDist(size_t i, size_t j) {return sqrt((getX(i) - getX(j))*(getX(i) - getX(j)) + (getY(i) - getY(j))*(getY(i) - getY(j)));};
        MPI_Win* getWin(unsigned winIndex) {return &win[winIndex];};
        void pushOffset(size_t globalNode);

        // INTERFACE
        void createNodes(setting* settings);
        void replaceTemp();
        void exchangeTemp();
        void accumulateTemp();
        void exchangeMass();
        void accumulateMass();
        void exchangeBCType();
        void accumulateBCType();

        // LOL
        void printBoundaryConds();
};

#endif /* NODE_H_ */
