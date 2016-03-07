#include "node.h"
#include "utility.h"
#include "mpi.h"
#include <vector>

/***************************************************************************************************
 * void node::createNodes()
 * This routine calls private methods to read mxyz and data file.
 **************************************************************************************************/
void node::createNodes(setting* settings)
{
    // Initialize the private variables with zero.
    for (int in = 0; in < nn; in++)
    {
        BCtype[in] = 0;
        x[in]      = 0.0;
        y[in]      = 0.0;
        T[in]      = 0.0;
        mass[in]   = 0.0;
    }
  
    readMxyzPar(settings);
    readDataPar(settings);
    myRank = settings->getRank();
    numProcs = settings -> getProcs();
    maxNodes = (settings->getNn() - 1) / numProcs + 1;
    /// Declare vectors here 
    
    // cout << "nodecpp:YOU HAVE NO IDEA HOW VECTORS WORK!\n";
    offsets.resize(numProcs+1);
    //localOffsets.resize(numProcs);
    // cout << "nodecpp:OK, YOU DO!\n";
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}


void node::initMPI()
{
    sortOffsets();
    assignRemoteStorage();
    createWindows();
    createDatatypes();
    exchangeData();
}

/***************************************************************************************************
 * void node::readMxyz()
 ***************************************************************************************************
 * READ THE MXYZ FILE
 * mxyz file contains the node coordinates
 **************************************************************************************************/
void node::readMxyz(setting* settings)
{
    int         nn;             // nuber of nodes
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files

    nn    = settings->getNn();
    dummy = settings->getMxyzFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    readStream = new char [nsd*sizeof(double)];
    file.seekg (0, ios::beg);
    for(int i=0; i<nn; i++)
    {
        file.read(readStream, nsd*sizeof(double));
        utility::swapBytes(readStream, nsd, sizeof(double));
        x[i] = *((double*)readStream);
        y[i] = *((double*)readStream+1);
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();

    return;
}

/*****
 * Parallel read for mxyz file
 * ****/
void node::readMxyzPar(setting* settings)
{
  MPI_File fd;
  MPI_Offset offset;

  if (settings -> getRank() == settings -> getProcs() - 1)
    offset = sizeof(double)*(settings -> getNn() - getNn())*nsd;
  else
    offset = sizeof(double)*settings -> getRank()* getNn()*nsd;
  
  MPI_File_open(MPI_COMM_WORLD, settings -> getMxyzFile().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);

  if (fd == NULL)
  {
    cout << "Unable to open file : " << settings -> getMxyzFile() << "\n";
    exit(0);
  }

  double* readTmp = new double [nsd*getNn()];
  MPI_File_set_view(fd, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

  MPI_File_read(fd, readTmp, nsd*getNn(), MPI_DOUBLE, MPI_STATUS_IGNORE);
  utility::swapBytes((char *)readTmp, nsd*getNn(), sizeof(double));
  
  for(int i = 0; i < getNn(); ++i)
  {
    x[i] = readTmp[2*i];
    y[i] = readTmp[2*i + 1];
  }
  
  delete[] readTmp;
  MPI_File_close(&fd);
//  printf("Proc #%d MXYZ done.\n", settings -> getRank());  

  // printf("proc #%d done with mxyz file\n Has %ld nodes!", settings -> getRank(), getNn()); 

  // cout << "I am proc: # " << settings -> getRank() << "\n";
  // for (int i = 0; i < getNn(); ++i)
  //  cout << x[i] << " " << y[i] << "\n";
 
  // exit(0); 
  
  return;
}

/***************************************************************************************************
 * void node::readData()
 ***************************************************************************************************
 * READ THE DATA FILE
 * data file contains the initial value of the field.
 **************************************************************************************************/
void node::readData(setting* settings)
{
    int      nn;             // number of nodes to be determined from minf file
    bool     restart;        // Is this a restart?
    ifstream file;           // file name obj
    string   dummy;          // dummy string to hold names
    char*    readStream;     // temperory var used for strings read from files
    double   dummyDouble;    // temperory var used for double values read from files

    nn      = settings->getNn();
    restart = settings->getRestart();
  
    if (restart)
    {
        dummy = settings->getDataFile();
        file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
        if (file.is_open()==false)
        {
            cout << "Unable to open file : " << dummy << endl;
            exit(0);
        }
        readStream = new char [sizeof(double)];
        file.seekg (0, ios::beg);
        for(int i=0; i<nn; i++)
        {
            file.read (readStream, sizeof(double));
            utility::swapBytes(readStream, 1, sizeof(double));
            T[i] = *((double*)readStream);
        }
        cout << "> File read is completed: " << dummy << endl;
        file.close();
    }
    else
    {
        dummyDouble = settings->getInitT();
        for(int i=0; i<nn; i++)
            T[i] = dummyDouble;
        cout << "> Data initialization is completed." << endl;
    }

    return;
}

void node::readDataPar(setting* settings)
{
    int      nn;             // number of nodes to be determined from minf file + adjusted to current number of Procs
    bool     restart;        // Is this a restart?
    ifstream file;           // file name obj
    string   dummy;          // dummy string to hold names
    char*    readStream;     // temperory var used for strings read from files
    double   dummyDouble;    // temperory var used for double values read from files

    nn      = getNn();       // LOCAL NUMBER !!!
    restart = settings->getRestart();

    /* Currently restart option is unavailable !! */
    if (restart) // Do nothing!
    {
    //         //     dummy = settings->getDataFile();
    //         //     file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    //         //     if (file.is_open()==false)
    //         //     {
    //         //         cout << "Unable to open file : " << dummy << endl;
    //         //         exit(0);
    //         //     }
    //         //     readStream = new char [sizeof(double)];
    //         //     file.seekg (0, ios::beg);
    //         //     for(int i=0; i<nn; i++)
    //         //     {
    //         //         file.read (readStream, sizeof(double));
    //         //         utility::swapBytes(readStream, 1, sizeof(double));
    //         //         T[i] = *((double*)readStream);
    //         //     }
    //         //     cout << "> File read is completed: " << dummy << endl;
    //         //     file.close();
    }
    else
    {
      dummyDouble = settings->getInitT();
      for(int i=0; i<nn; i++)
        T[i] = dummyDouble;
 //     cout << "> Data initialization is completed." << endl;
    }
    return;
}

/***************************************************************************************************
 * void node::printBoundaryConds()
 * This routine prints all of those pesky boundary conditions to std out.
 **************************************************************************************************/
void node::printBoundaryConds()
{
    for(size_t in = 0; in < nn; ++in)
	std::cout << in << "\t" <<  BCtype[in] << std::endl; 

    return;
}

/***************************************************************************************************
 * void node::pushOffset()
 * This routine stores the positions of nodal values on the remote processes.
 **************************************************************************************************/
void node::pushOffset(size_t globalNode)
{
     int process = getRemoteRank(globalNode);
     
     if (process == myRank) return;
     

     size_t node = getRemoteNode(globalNode);
     offsets[process].push_back(node);
     //localOffsets[process].push_back(additionalNodes);

     return;
}

/***************************************************************************************************
* void node::sortOffset()
* This routine sorts and removes duplicates from the offsets vector.
***************************************************************************************************/
void node::sortOffsets()
{
    for (int i = 0; i < numProcs; ++i)
    {
        std::sort(offsets[i].begin(), offsets[i].end());
        offsets[i].erase( std::unique(offsets[i].begin(), offsets[i].end()), offsets[i].end() );
        offsets[numProcs].push_back(offsets[i].size());
        additionalNodes += offsets[i].size();
    }

    for (int i = numProcs-1; i > 0; --i)
        offsets[numProcs][i] = offsets[numProcs][i-1];

    offsets[numProcs][0] = 0;

    for (int i = 1; i < numProcs; ++i)
        offsets[numProcs][i] += offsets[numProcs][i-1];  

    return;
}    
/***************************************************************************************************
* void node::assignRemoteStorage()
* This routine assigns space for the local storage of remote variables
***************************************************************************************************/
void node::assignRemoteStorage()
{
   // utility::resizeVec(remoteBC, offsets);
   // utility::resizeVec(remoteX, offsets);
   // utility::resizeVec(remoteY, offsets);
   // utility::resizeVec(remoteT, offsets);
   // utility::resizeVec(remoteM, offsets);
    BCtype.resize(nn + additionalNodes);
    x.resize(nn + additionalNodes);
    y.resize(nn + additionalNodes);
    T.resize(nn + additionalNodes);
    mass.resize(nn + additionalNodes);

    return;
}

/***************************************************************************************************
* void node::createDatatypes()
* This routine creates the derived datatypes used for communication
***************************************************************************************************/
void node::createDatatypes()
{
    remoteDatatypeInt = new MPI_Datatype[numProcs];
    remoteDatatype = new MPI_Datatype[numProcs];
    localDatatypeInt = new MPI_Datatype[numProcs];
    localDatatype = new MPI_Datatype[numProcs];
    for(int i = 0; i < numProcs; ++i)
    {
//        if (i == myRank) continue;        

        MPI_Type_create_indexed_block(offsets[i].size(), 1, &offsets[i][0], MPI_INT, &remoteDatatypeInt[i]);
        MPI_Type_create_indexed_block(offsets[i].size(), 1, &offsets[i][0], MPI_DOUBLE, &remoteDatatype[i]);
        MPI_Type_contiguous(offsets[i].size(), MPI_INT, &localDatatypeInt[i]);
        MPI_Type_contiguous(offsets[i].size(), MPI_DOUBLE, &localDatatype[i]);
    
        MPI_Type_commit (&remoteDatatypeInt[i]);
        MPI_Type_commit (&remoteDatatype[i]);
        MPI_Type_commit (&localDatatypeInt[i]);
        MPI_Type_commit (&localDatatype[i]);
    }
    return;
}
/***************************************************************************************************
* void node::deleteDatatypes()
* This routine frees the derived datatypes used for communication
***************************************************************************************************/
void node::deleteDatatypes()
{
    delete[] remoteDatatypeInt;
    delete[] remoteDatatype;
    delete[] localDatatypeInt;
    delete[] localDatatype;
    return;
}

/***************************************************************************************************
* void node::createWindows()
* This routine creates the memory windows required for RMA communication
***************************************************************************************************/
void node::createWindows()
{
    MPI_Win_create(&BCtype[0], nn*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win[0]);
    MPI_Win_create(&x[0], nn*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win[1]);
    MPI_Win_create(&y[0], nn*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win[2]);
    MPI_Win_create(&T[0], nn*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win[3]);
    MPI_Win_create(&mass[0], nn*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win[4]);
    return;
}


/***************************************************************************************************
* void node::freeWindows()
* This routine frees the memory windows used for communication
***************************************************************************************************/
void node::freeWindows(int numWindows)
{
    for (int i = 0; i < numWindows; ++i)  MPI_Win_free(&win[i]);
    return;
}

/***************************************************************************************************
* void node:exchangeData()
* This routine exchanges the initial nodal data required to construct the local matrices EXCLUDES TEMP 
***************************************************************************************************/
void node::exchangeData()
{
//    for (int i = 0; i < numProcs; ++i)
//    {
        // Skip own rank, no attempt to 'get' from self
//        if (i == myRank) continue;
        // Get Boundary Conditions
/*        MPI_Win_fence(0, win[0]);
        MPI_Get(&BCtype[nn + getExtra(i)], 1, localDatatypeInt[i], 
                       i, 0, 1, remoteDatatypeInt[i], win[0]);
        MPI_Win_fence(0, win[0]);
*/        // Get x-coordinates
        MPI_Win_fence(0, win[1]);
        for (int i = 0; i < numProcs; ++i)
            MPI_Get(&x[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[1]);
        MPI_Win_fence(0, win[1]);
        // Get y-coordniates
        MPI_Win_fence(0, win[2]);
        for (int i = 0; i < numProcs; ++i)
            MPI_Get(&y[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[2]);
        MPI_Win_fence(0, win[2]);
/*        // Get Temperature
        MPI_Win_fence(0, win[3]);
        MPI_Get(&T[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[3]);
        MPI_Win_fence(0, win[3]);
        // Get Mass
        MPI_Win_fence(0, win[4]);
        MPI_Get(&mass[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[4]);
        MPI_Win_fence(0, win[4]);
*/
//    }
    return;
}



/***************************************************************************************************
* void node:replaceTemp()
* This routine replaces the local nodal temperature with values from remote procs 
***************************************************************************************************/
void node::replaceTemp()
{
    MPI_Win_fence(0, win[3]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
        // if (i == myRank) continue;
        // Get Temperature
//        MPI_Win_fence(0, win[3]);
        if (i != myRank)
            MPI_Accumulate(&T[nn + getExtra(i)], 1, localDatatype[i],
                              i, 0, 1, remoteDatatype[i], MPI_REPLACE, win[3]);
//        MPI_Win_fence(0, win[3]);
    }
    MPI_Win_fence(0, win[3]);
    return;
}
/***************************************************************************************************
* void node:exchangeTemp()
* This routine exchanges the nodal temperature data required at every time step 
***************************************************************************************************/
void node::exchangeTemp()
{
    MPI_Win_fence(0, win[3]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
 //       if (i == myRank) continue;
        // Get Temperature
//        MPI_Win_fence(0, win[3]);
        MPI_Get(&T[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[3]);
//        MPI_Win_fence(0, win[3]);
    }
    MPI_Win_fence(0, win[3]);
    return;
}

/***************************************************************************************************
* void node:accumulateTemp()
* This routine acumulates the temperatures calculated on remote nodes to it's local storage 
***************************************************************************************************/
void node::accumulateTemp()
{
    MPI_Win_fence(0, win[3]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
 //       if (i == myRank) continue;
        // Accumulate Temperature
//        MPI_Win_fence(0, win[3]);
        MPI_Accumulate(&T[nn + getExtra(i)], 1, localDatatype[i],
                              i, 0, 1, remoteDatatype[i], MPI_SUM, win[3]);
//        MPI_Win_fence(0, win[3]);
    }
    MPI_Win_fence(0, win[3]);
    return;
}




/***************************************************************************************************
* void node:exchangeMass()
* This routine exchanges the nodal mass 
***************************************************************************************************/
void node::exchangeMass()
{
    MPI_Win_fence(0, win[4]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
 //       if (i == myRank) continue;
        // Get Temperature
//        MPI_Win_fence(0, win[4]);
        MPI_Get(&mass[nn + getExtra(i)], 1, localDatatype[i], 
                       i, 0, 1, remoteDatatype[i], win[4]);
//        MPI_Win_fence(0, win[4]);
    }
    MPI_Win_fence(0, win[4]);
    return;
}

/***************************************************************************************************
 * * void node:accumulateMass()
 * * This routine acumulates the mass calculated on remote nodes to it's local storage 
 * ***************************************************************************************************/
void node::accumulateMass()
{ 
    MPI_Win_fence(0, win[4]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
        // if (i == myRank) continue;
        // Accumulate Mass
//        MPI_Win_fence(0, win[4]);
        MPI_Accumulate(&mass[nn + getExtra(i)], 1, localDatatype[i],
                              i, 0, 1, remoteDatatype[i], MPI_SUM, win[4]);
//        MPI_Win_fence(0, win[4]);
    }
    MPI_Win_fence(0, win[4]);
    return;
}


/***************************************************************************************************
* void node:exchangeBCType()
* This routine exchanges the boundary conditions 
***************************************************************************************************/
void node::exchangeBCType()
{
    MPI_Win_fence(0, win[0]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
 //       if (i == myRank) continue;
        // Get Temperature
//        MPI_Win_fence(0, win[0]);
        MPI_Get(&BCtype[nn + getExtra(i)], 1, localDatatypeInt[i], 
                       i, 0, 1, remoteDatatypeInt[i], win[0]);
//        MPI_Win_fence(0, win[0]);
    }
    MPI_Win_fence(0, win[0]);
    return;
}
/***************************************************************************************************
 * * void node:accumulateBCType()
 * * This routine places the max boundary condition value into local storage 
 * ***************************************************************************************************/
void node::accumulateBCType()
{
    MPI_Win_fence(0, win[0]);
    for (int i = 0; i < numProcs; ++i)
    {
        // Skip own rank, no attempt to 'get' from self
        // if (i == myRank) continue;
        // Accumulate Mass
//        MPI_Win_fence(0, win[0]);
        MPI_Accumulate(&BCtype[nn + getExtra(i)], 1, localDatatypeInt[i],
                              i, 0, 1, remoteDatatypeInt[i], MPI_MAX, win[0]);
//        MPI_Win_fence(0, win[0]);
    }
    MPI_Win_fence(0, win[0]);
    return;
}

/*
 * Returns index Should be done for the offset array!!! 
 * */

size_t node::searchArr(size_t in)
{
    size_t remoteNode = getRemoteNode(in);
    int remoteRank = getRemoteRank(in);
    for (int i =0; i < offsets[remoteRank].size(); ++i)
    {
        if (remoteNode == offsets[remoteRank][i])
        return i;
    }
}

