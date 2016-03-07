#include "element.h"
#include "utility.h"
#include "mpi.h"

/***************************************************************************************************
 * void element::createElements()
 * This routine calls private methods to read mien and mrng file.
 **************************************************************************************************/
void element::createElements(setting* settings, node* nodes)
{
    // Initialize all private variables
    for (int i = 0; i < ne*nen; i++)
    {
        conn[i] = 0;
        M[i]    = 0.0;
        F[i]    = 0.0;
        RHS[i] = 0.0;
        bound[i] = 0.0;
    }
    for (int i = 0; i < ne*nef;     i++)    FG[i]  = 0;
    for (int i = 0; i < ne*nen*nen; i++)    K[i]   = 0.0; 

    maxElems = (settings->getNe() - 1) / settings->getProcs() + 1;
    // printf("%d: init finished createElements(Before)\n", settings -> getRank());
    // Start element operations
    readMienPar(settings);
    // readMien(settings);
    readMrngPar(settings);
    
    calculateOffsets(nodes);
//    nodes -> sortOffsets();

    MPI_Barrier(MPI_COMM_WORLD); // KEEP THIS ONE 
    nodes -> initMPI();    
    modifyConn(nodes);
    calculateJacobian(nodes);
    calculateElementMatrices(settings, nodes);
    nodes->accumulateMass();
    nodes->exchangeMass();
    applyBoundaryConditions(settings, nodes);

    return;
}


/***************************************************************************************************
 * void element::readMien()
 * This routine reads mien file.
 **************************************************************************************************/
void element::readMien(setting* settings)
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files

    dummy = settings->getMienFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [nen*sizeof(int)];
    file.seekg (0, ios::beg);
    for(size_t i=0; i<ne; i++)
    {
        file.read (readStream, nen*sizeof(int));
        utility::swapBytes(readStream, nen, sizeof(int));
        for(int j=0; j<nen; j++)
            conn[i*nen+j] = *((int*)readStream+j)-1;
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();

    return;
}

/*****
 parallel read for mien file
****/ 
void element::readMienPar(setting* settings)
{
  MPI_File fd;
  MPI_Offset offset;

/* 
  ostringstream int2str; 
  int2str << setfill('0') << setw(5) << settings -> getRank();
  string strong = int2str.str();
  string aling;
  if (settings -> getMprmFile().empty())
  {
      // do nothing ;
  }
  else 
  {
      aling = settings -> getMprmFile().append(strong);
 //   int2str << getMprmFile().c_str() << setfill('0') << setw(5) << settings -> getRank();
  }   
*/
  if (settings -> getRank() == settings -> getProcs() - 1)
    offset = sizeof(int)*(settings -> getNe() - getNe()) *nen;
  else
    offset = sizeof(int)*settings -> getRank()* getNe()*nen;
//  int* permutation = new int[getNe()];
//  MPI_File_open(MPI_COMM_WORLD, aling.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);
//  MPI_File_set_view(fd, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//  MPI_File_read(fd, permutation, getNe(), MPI_INT, MPI_STATUS_IGNORE);
//  MPI_File_close(&fd);
//  utility::swapBytes((char *)permutation, ne, sizeof(int));
//  for (size_t ie = 0; ie < ne; ++ie)
//      --permutation[ie];

  MPI_File_open(MPI_COMM_WORLD, settings -> getMienFile().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);

  if (fd == NULL)
  {
    cout << "Unable to open file : " << settings -> getMienFile() << "\n";
    exit(0);
  }

  MPI_File_set_view(fd, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

  MPI_File_read(fd, conn, nen*getNe(), MPI_INT, MPI_STATUS_IGNORE);
/*
  for (size_t i = 0; i < getNe(); ++i)
  {
      offset = sizeof(int)*(permutation[i]*nen);
      MPI_File_set_view(fd, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
      MPI_File_read(fd, &conn[i*nen], nen, MPI_INT, MPI_STATUS_IGNORE);
  }

  delete [] permutation;
*/
  utility::swapBytes((char *)conn, nen*getNe(), sizeof(int));
  for (size_t ie = 0; ie < ne; ++ie)
    for (int ien = 0; ien < nen; ++ien)
    {
      --conn[ie*nen+ien];
    }

  MPI_File_close(&fd);

//  printf("Proc #%d MIEN done.\n", settings -> getRank()); 

  // cout << "I am proc: # " << settings -> getRank() << "\n";
  // cout << "those should be different!" << getNe() << " " << settings -> getNe() << "\n";
  //for (int i = 0; i < getNe()*nen; i += nen)
  //  cout << conn[i] << " " << conn[i + 1] << " " <<  conn[i + 2] << "\n";

  //exit(0);

  return;
}

/***************************************************************************************************
 * void element::readMrng()
 * This routine reads mrng file.
 **************************************************************************************************/
void element::readMrng(setting* settings)
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files

    dummy = settings->getMrngFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [nef*sizeof(int)];
    file.seekg (0, ios::beg);
    for(size_t i=0; i<ne; i++)
    {
        file.read (readStream, nef*sizeof(int));
        utility::swapBytes(readStream, nef, sizeof(int));
        for(ptrdiff_t j=0; j<nef; j++)
            FG[i*nef+j] = *((int*)readStream+j);
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();

    return;
}

void element::readMrngPar(setting* settings)
{
  MPI_File fd;
  MPI_Offset offset;

  if (settings -> getRank() == settings -> getProcs() - 1)
    offset = sizeof(int)*(settings -> getNe() - getNe())*nef;
  else
    offset = sizeof(int)*settings -> getRank()* getNe()*nef;

  MPI_File_open(MPI_COMM_WORLD, settings -> getMrngFile().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);

  if (fd == NULL)
  {
    cout << "Unable to open file : " << settings -> getMrngFile() << "\n";
    exit(0);
  }

  MPI_File_set_view(fd, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

  MPI_File_read(fd, FG, nef*getNe(), MPI_INT, MPI_STATUS_IGNORE);
  utility::swapBytes((char *)FG, nef*getNe(), sizeof(int));

  MPI_File_close(&fd);

//  printf("Proc #%d MRNG done.\n", settings -> getRank());  
  /* cout << "I am proc: # " << settings -> getRank() << "\n";
  //   // cout << "those should be different!" << getNe() << " " << settings -> getNe() << "\n";
  for (int i = 0; i < getNe()*nef; i += nef)
      cout << FG[i] << " " << FG[i + 1] << " " <<  FG[i + 2] << "\n";
  
  exit(0);
  */
  return;
  }


/***************************************************************************************************
 * void element::calculateJacobian()
 * This routine calculates the Jacobian determinant and derivatives for each GQP.
 **************************************************************************************************/
void element::calculateJacobian(node* nodes)
{
    for (int ie = 0; ie < ne; ++ie){ // for each element 

        // if (ie == 0){ //check point
        //     cout << "THIS LINE IS IMPORTANT\n" << nodes -> getX(getConn(0, 0)) << " " << nodes -> getX(getConn(0, 1)) << " " << nodes -> getX(getConn(0, 2)) << "\n";
        //     cout << "THIS LINE IS IMPORTANT\n" << nodes -> getY(getConn(0, 0)) << " " << nodes -> getY(getConn(0, 1)) << " " << nodes -> getY(getConn(0, 2)) << "\n";
        // }

        double dummyJ[4]; // dummy Jacobi matrix 11 12 21 22  
        double dummy = 0; 
        for (int iGQP = 0; iGQP < getNGQP(); ++iGQP){ //for each quadrature point
            dummyJ[0] = dummyJ[1] = dummyJ[2] = dummyJ[3] = 0;
            for (int ien = 0; ien < getNen(); ++ien){
                dummyJ[0] += getDSdKsi(iGQP, ien) * nodes -> getX(getConn(ie, ien));
                dummyJ[1] += getDSdKsi(iGQP, ien) * nodes -> getY(getConn(ie, ien));
                dummyJ[2] += getDSdEta(iGQP, ien) * nodes -> getX(getConn(ie, ien));
                dummyJ[3] += getDSdEta(iGQP, ien) * nodes -> getY(getConn(ie, ien));
            }
            dummy = dummyJ[0]*dummyJ[3] - dummyJ[2]*dummyJ[1];
            setDetJ(ie, iGQP, dummy); // set the Jacobian 

            // cout << "\n JACOBIAN = " << getDetJ(ie, iGQP) <<"\n";

            /* Find an inverse of Jacobi matrix */
	    //cout << "Before: " << dummyJ[0] << " " << dummyJ[3] << " ";
            
	    swap (dummyJ[0], dummyJ[3]); // swap elements on the main diagonal
            //cout << "After: " << dummyJ[0] << " " << dummyJ[3] << "\n";
	    dummyJ[0] /= getDetJ(ie, iGQP);
            dummyJ[1] /= -getDetJ(ie, iGQP);
            dummyJ[2] /= -getDetJ(ie, iGQP);
            dummyJ[3] /= getDetJ(ie, iGQP);
            /* Now dummyJ stores the inverse Jacobi matrix */
            
            /* Assign the same values to each DSdx in a row, 
                    propably WRONG*/
            for (int ien = 0; ien < getNen(); ++ien){
                setDSdx(ie, iGQP, ien, getDSdKsi(iGQP, ien)*dummyJ[0] + getDSdEta(iGQP, ien)*dummyJ[1]);
                setDSdy(ie, iGQP, ien, getDSdKsi(iGQP, ien)*dummyJ[2] + getDSdEta(iGQP, ien)*dummyJ[3]);       
            }      
        }
    }

    return;
}


/***************************************************************************************************
 * void element::calculateElementMatrices()
 * This routine calculates element mass, stiffness and forcing matrices
 * as well as nodal mass.
 **************************************************************************************************/
void element::calculateElementMatrices(setting* settings, node* nodes)
{
    size_t totNodes = nodes->getNn() + nodes->getAdditionalNodes();
    //Setting the nodal mass to zero
    for (int in = 0; in < totNodes; ++in)
        nodes->setMass(in, 0);

    for (size_t ie = 0; ie < ne; ++ie)
    {
    	double* consistentM = new double[nen*nen];
        for (int ien = 0; ien < nen; ++ien)
        {
            setForcingElement(ie, ien, settings->getSource());
            for (int jen = 0; jen < nen; ++jen)
            {
                setStiffnessElement(ie, ien, jen, settings->getD());     
                setConsistentMassElement(ie, ien, jen, &consistentM[0]);
            }
        }
        setLumpedMassMatrix(ie, &consistentM[0]);
        delete[] consistentM;
//        printElementStiffnessMatrix(ie);
//        printElementMassMatrix(ie);
     
         // Summing the element contributions to each node mass   
         for (int ien = 0; ien < nen; ++ien)
         {
     	    nodes->addMass(getConn(ie, ien), M[ie*nen+ien]); 
         }  
     }
//for (size_t i = 0; i < nodes->getNn(); ++i) cout<<nodes->getMass(i)<<endl;

    return;
}


/***************************************************************************************************
 * void element::applyBoundaryConditions()
 * Applies Drichlet and Neumann boundary conditions.
 **************************************************************************************************/
void element::applyBoundaryConditions(setting* settings, node* nodes)
{

// PROBABLY WE NEED SOME SEARCH FOR THE CORRECT FACE ;)

    // RHS assembly + BC  
    for (size_t ie = 0; ie < ne; ++ie)
    {
        for (size_t ien = 0; ien < nen; ++ien)
        {
	        size_t localNode = ie*nen + ien;
                size_t globalNode = getConn(ie, ien);
                int faceGroup = getFG(ie, ien);
/* Here we set boundary conditions based on the global face groups. We skip internal faces denoted*/
/* by faceGroup = 0										  */
                if (faceGroup > 0)
                { 
/* Von Neumann and Robin boundaries are additive, ie each element provides a portion (determined  */
/* by the shape funtions) of the boundary condition to each of its nodes. However Dirichelet      */
/* conditions are not additive, the boundary value is assigned rather than incremented. This      */
/* complicates the assignment of boundary values on the nodes where Dirichelet and Neumann (or    */
/* Robin) faces meet. The two if statements below prevent Dirichelet conditions from overwriting  */
/* already assigned Neumann and Robin boundary values                                             */ 
                    int nodeFace = abs(nodes->getBCtype(globalNode));
                    int currentType = settings->getBC(nodeFace)->getType();
                    int proposedType = settings->getBC(faceGroup)->getType();
                    //std::cout << settings->getBC(nodes->getBCtype(4))->getType() << std::endl;
                    if ((currentType > 1) & (proposedType == 1))
                    {
                        continue;
                    }
/* Set the nodes either side of the current face to the face group of the current face. Positive */
/* value indicates the node and face have the same local index. Negative sign denotes the node   */
/* has the next local index ie (ien+1)%nen. This is required to calculate the length of face.     */
                    nodes->setBCtype(globalNode, faceGroup);
                    if (proposedType > 1)
                    {
                        nodes->setBCtype(getConn(ie, (ien+1)%(nen)), faceGroup);
//std::cout<<getConn(ie, (ien+1)%nen) << "\t" <<settings->getBC(abs(nodes->getBCtype(getConn(ie,(ien+1)%nen))))->getType() << std::endl;;
                    }

                    
                    if (proposedType == 2)
                    {
                        double distance = nodes->getDist(globalNode, getConn(ie, (ien+1)%nen));
                        setBound(ie, ien, distance);
                        setBound(ie, (ien+1)%nen, -distance); 
                    }
                    else if (proposedType == 3)
                    {
                        double distance = nodes->getDist(globalNode, getConn(ie, (ien+1)%nen));
                        setBound(ie, ien, distance);
                        setBound(ie, (ien+1)%nen, -distance);
                    }




       
                } 
        }            
    }

    nodes->accumulateBCType();
    nodes->exchangeBCType();

    for (size_t ie = 0; ie < ne; ++ie)
    {
         for (size_t ien = 0; ien < nen; ++ien)
         {
         // RESET THE K STIFFNESS VALUES
         // Put everything to zero except elements on the main diagonal
              size_t globalNode = getConn(ie,ien);
              size_t localNode = ie*nen + ien;
/* Extract the face group encoded globally for the local node. Determine its sign. Use the       */
/* absolute value of the face group (it is an index into an array)                               */
              int nodeFace = nodes->getBCtype(globalNode);
              int boundaryType = settings->getBC(nodeFace)->getType();


/* If the face is a Dirichelet face, find the value stored in the settings object assign it to   */
/* the node and modify the elementwise stiffness matrix                                          */
              if (boundaryType == 1) 
              {
                    double dirichlet = settings->getBC(nodeFace)->getValue1(); 
                    setBound(ie, ien, dirichlet);
                    nodes->setT(globalNode, dirichlet);
                    for (size_t jen = 0; jen < nen; ++jen)
                    { // CAN BE EROUNEOUS
                        if (jen != ien)
                            setK(ie, ien, jen, 0.0);
                        else
                            setK(ie, ien, jen, 1.0);
                    }
              }
/* If the face is a Neumann face, find the flux value stored in settings and set the boundary val */
/* to the heat flow across half the face                                                          */
              else if (boundaryType == 2)
              {
/* Calculate the length of the face. Here the value 'dir' encodes the direction to the node on   */
/* other end of the face                                                                         */
                  double distance = getBound(ie, ien);//nodes->getDist(globalNode, getConn(ie, (ien+dir)%nen));
                  int dir = (distance > 0) - (distance < 0);
                  if (0 == dir) continue;
                  distance = abs(distance);
                  double halfHeat = settings->getBC(nodeFace)->getValue1() * 0.5;
                  setBound(ie, ien, distance * halfHeat);

              }
/* If the face is a Robin face, find the temperature of the resevoir and the transfer coefficent */
/* Set the boundary value and modify the elementwise stiffness matrix                            */
              else if (boundaryType == 3)
              {         
/* Calculate the length of the face. Here the value 'dir' encodes the direction to the node on   */
/* other end of the face                                                                         */
                  double distance = getBound(ie, ien);//nodes->getDist(globalNode, getConn(ie, (ien+dir)%nen));
                  int dir = (distance > 0) - (distance < 0);
                  if (0 == dir) continue;
                  distance = abs(distance);
                  double transferCoeff = settings->getBC(nodeFace)->getValue2();
                  double temperature = settings->getBC(nodeFace)->getValue1();   
                  setBound(ie, ien, distance * transferCoeff * temperature * 0.5);


                  for (size_t jen = 0; jen < nen; ++jen)
                  {
                      if (jen == ien)
                          addK(ie, ien, jen, transferCoeff * distance/3);
                      else if (jen == (ien + dir)%nen)
                          addK(ie, ien, jen, transferCoeff * distance/6);
                  }
              }
          }
 //     std::cout << "Element " << ie << std::endl;
 //     printElementStiffnessMatrix(ie);
      }
      nodes->replaceTemp();
//      nodes->exchangeTemp();
      return;
}


/***************************************************************************************************
 * void element::setStiffnessElement()
 * Uses Gauss quadrature to calculate one value of the element wise stiffness matrix.
 **************************************************************************************************/
void element::setStiffnessElement(size_t ie, size_t ien, size_t jen, double diffuse)
{

    double stiffnessElement = 0;
    
    for (int iGQP = 0; iGQP < nGQP; ++iGQP) 
    {
        double DSdx_i = getDSdx(ie, iGQP, ien);
        double DSdx_j = getDSdx(ie, iGQP, jen);
        double DSdy_i = getDSdy(ie, iGQP, ien);
        double DSdy_j = getDSdy(ie, iGQP, jen);

        double weight = getWeight(iGQP);
        double detJ = getDetJ(ie, iGQP);
        
        stiffnessElement += (DSdx_i*DSdx_j + DSdy_i*DSdy_j)*weight*detJ;
    }
    setK(ie, ien, jen, diffuse*stiffnessElement);

    return;
}

/***************************************************************************************************
 * void element::setForcingElement()
 * Uses the source term in setting.h to calculate a meaningless forcing element
 **************************************************************************************************/
void element::setForcingElement(size_t ie, size_t ien, double f)
{
    double forcingElement = 0;
    for (int iGQP = 0; iGQP < nGQP; ++iGQP) 
    {
        double S_i = getS(iGQP, ien);

        double weight = getWeight(iGQP);
        double detJ = getDetJ(ie, iGQP);
        
        forcingElement += S_i*weight*detJ;
    }
    setF(ie, ien, f*forcingElement);
    
    return;
}

/***************************************************************************************************
 * void element::setConsistentMassElement()
 * Uses Gauss quadrature to calculate one value of the element wise consistent mass matrix.
 **************************************************************************************************/
void element::setConsistentMassElement(size_t ie, size_t ien, size_t jen, double* consistentMassMatrix)
{

    double consistentMassElement = 0;
    
    for (int iGQP = 0; iGQP < nGQP; ++iGQP)
    {
      	double S_i = getS(iGQP, ien);
       	double S_j = getS(iGQP, jen);
    	double weight = getWeight(iGQP);
    	double detJ = getDetJ(ie, iGQP);
       	consistentMassElement += S_i * S_j * weight * detJ;
    }
    consistentMassMatrix[ien*nen + jen] = consistentMassElement;
    return;
}

/***************************************************************************************************
 * void element::setLumpedMassMatrix()
 * Lumps the element mass matrix 'consistentMassMatrix' 
 **************************************************************************************************/
void element::setLumpedMassMatrix(size_t ie, double* consistentMassMatrix)
{
    
    double massOfElement = 0;
    double massOfDiagonal = 0;
    double* lumpedMassMatrix = new double[nen];
    for (int ien = 0; ien < nen; ++ien)
    {
    	massOfDiagonal += consistentMassMatrix[ien*nen + ien];
    	for (int jen = 0; jen < nen; ++jen)
    	{
    		massOfElement += consistentMassMatrix[ien*nen + jen];
    	}
    	lumpedMassMatrix[ien] = consistentMassMatrix[ien*nen + ien];
    }
    
    for (int ien = 0; ien < nen; ++ien)
    { 
    	lumpedMassMatrix[ien] *= (massOfElement / massOfDiagonal);
    	setM(ie, ien, lumpedMassMatrix[ien]);
    }   
    
    delete[] lumpedMassMatrix;
    return;
}

/***************************************************************************************************
 * void element::printElementStiffnessMatrix()
 * Outputs the stiffness matrix to stdout
 **************************************************************************************************/
void element::printElementStiffnessMatrix(size_t ie)
{
	for (int ien = 0; ien < nen; ++ien)
    {
    	for (int jen = 0; jen < nen; ++jen)
    	{
    		std::cout << K[ie*nen*nen+ien*nen+jen] << "\t";
    	}
    	std::cout << std::endl;
    }
    std::cout << std::endl;
    return;
}

// void element::SomeRandomFunctionToCollectRHS(/*with some random args*/)
// {
//     double localRHS[nen];
//     memset(localRHS, 0, nen*sizeof(double));

//     for (int ie = 0; ie < ne; ++ie){
//         memset(localRHS, 0, nen*sizeof(double));
//         for (int ien; ien < nen; ++ien){
//             localRHS[ie] = dt*(F[ien*ie] + /* Skip BC here */;
//         }
//         // flush localRHS to globalRHS[nn] matrix
//     }


//     return;
// }



/***************************************************************************************************
 *  * void element::calculateRHS()
 *  * Recalculates the RHS vector of the time dependent problem
 *  * This is called every time step
****************************************************************************************************/
void element::calculateRHS(node* nodes, setting* settings)
{
    //applyBoundaryConditions(settings, nodes);
    double dT = settings->getDt();
    double* RHS = getRHSptr(0);
    double* bound = getBoundptr(0);
    
    for(size_t ie = 0; ie < ne; ++ie)
    {
	for(int ien = 0; ien < nen; ++ien)
	{
	    size_t localNode = ie*nen + ien;
	    size_t globalNode =getConn(ie, ien);
	    RHS[localNode] = 0;
            RHS[localNode] += bound[localNode];
	    for(int jen = 0; jen < nen; ++jen)
	    {
	        double* K = getKptr(ie);
	        RHS[localNode] -=  K[ien*nen + jen] * nodes->getT(getConn(ie, jen));

	    }
	    double* F = getFptr(ie);
            RHS[localNode] += F[ien];
            RHS[localNode] *= dT;
            double* M = getMptr(ie);
            RHS[localNode] += M[ien] * nodes->getT(globalNode);
//            std::cout << localNode << "\t" << globalNode << "\t" << M[ien] << std::endl;            
	}
    }

    return;
}

/***************************************************************************************************
 * void element::printElementMassMatrix()
 * Outputs the mass matrix to stdout
 **************************************************************************************************/
void element::printElementMassMatrix(size_t ie)
{
    for (int ien = 0; ien < nen; ++ien)
    {
    	std::cout << ie*nen+ien << "\t" << getConn(ie, ien) << "\t" << M[ie*nen+ien] << "\n";
    }
    std::cout << std::endl;
    return;
}


/***************************************************************************************************
 * void element::calculateOffsets()
 * Calculates the offsets of the data in the remote processes
 **************************************************************************************************/
void element::calculateOffsets(node * nodes)
{
    for (size_t ie = 0; ie < ne; ++ie)
        for (int ien = 0; ien < nen; ++ien)
        {
            size_t globalNode = getConn(ie, ien);
            int rank = nodes->getRemoteRank(globalNode);
            int myRank = nodes->getRank();            
 //           if (rank == myRank)
 //               connCopy[ie*nen+ien] = nodes->getRemoteNode(globalNode);
            if (rank != myRank)
            {
                nodes->pushOffset(globalNode);
            }
        }       
    return;
}

/***************************************************************************************************
 * void element::modifyConn()
 * Changes connectivity array to account for remote nodes (which are appended to the end)
 **************************************************************************************************/
void element::modifyConn(node * nodes)
{
    int* tempConn = new int[ne*nen]();
    for (size_t ie = 0; ie < ne; ++ie)
        for (int ien = 0; ien < nen; ++ien)
        {
            size_t globalNode = getConn(ie, ien);
            size_t remoteNode = nodes->getRemoteNode(globalNode);
            int rank = nodes->getRemoteRank(globalNode);
            int myRank = nodes->getRank();
      //// SEARCH FOR AND RETURN index IN OFFSETS
            if (rank == myRank)
            {
                tempConn[ie*nen+ien] = remoteNode;
            }
            else if (rank != myRank)
            {
                size_t nn = nodes->getNn();
                int extra = nodes->getExtra(rank);
                int x = nodes->searchArr(globalNode);
                tempConn[ie*nen+ien] = nn + extra + x; 
            }
        }
    for (size_t i = 0; i < ne*nen; ++i) conn[i] = tempConn[i];
    delete[] tempConn;
    return;
}

