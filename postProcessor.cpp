#include "postProcessor.h"

/***************************************************************************************************
preProcessorControl
***************************************************************************************************/

void postProcessor::postProcessorControl(setting* argSettings, node* argNodes, element* argElements)
{
    int myRank, numProcs;             // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
 
    settings = argSettings;
    nodes = argNodes;
    elements = argElements;
 
//    if (myRank == 0) cout << endl << "================ POST-PROCESSING =================" << endl;

    //evaluateLimits();
    vtkVisualization();
    
    return;
}

/***************************************************************************************************
Evaluates the maximum and minimum temperatures in the field
***************************************************************************************************/
void postProcessor::evaluateLimits()
{
    int nn = settings -> getNn();
    int nnc = nodes -> getNn();
    int myRank, numProcs;             // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    double T;
    minT = std::numeric_limits<double>::max();
    maxT = std::numeric_limits<double>::min(); 

    for(int i = 0; i < nnc; ++i)
    {
      T = nodes -> getT(i + nodes->getMaxNodes()*myRank);
        if(T < minT)
            minT = T;
        if(T > maxT)
            maxT = T;
    }
//  cout << "Tmin " << minT << endl << "Tmax " << maxT << endl;

    return;
}

/***************************************************************************************************
// Main visualization function
***************************************************************************************************/
void postProcessor::vtkVisualization()
{
    static int count = 0;
    int nn = settings->getNn();
    int nnc = nodes->getNn();

//    int nnl = mesh->getNnl();

    int nna = nodes->getAdditionalNodes();
    int mnc = nodes->getMaxNodes();
    int ne = settings->getNe();
    int nec = elements->getNe();
    int mec = elements->getMaxElems();   
    int nen = elements->getNen();
    int myRank, numProcs;             // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    ostringstream int2str;
    int2str << myRank;
    string strMype = int2str.str(); 
    string dummy;
    int2str.str("");
    int2str << count;
    string tStep = int2str.str();

    // VTK Double Array
    vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
    pcoords->SetNumberOfComponents(3);
    pcoords->SetNumberOfTuples(nnc+nna);

    //vtkDoubleArray type pcoords is filled with the data in meshPoints.
    for (int i=0; i<nnc+nna; i++){
        pcoords -> SetTuple3(i, nodes->getX(i), nodes->getY(i), 0.0);
//        printf("Rank %d:\t i %d: %f %f\n", myRank, i + nodes -> getMaxNodes()*myRank, nodes->getX(i + nodes -> getMaxNodes()*myRank), nodes->getY(i + nodes -> getMaxNodes()*myRank));
    }
    

    //vtkPoints type outputPoints is filled with the data in pcoords.
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetData(pcoords);

    //Connectivity is written to vtkCellArray type outputCells
    vtkSmartPointer<vtkCellArray> connectivity = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i<nec; i++)
    {
        connectivity->InsertNextCell(nen);
        for(int j=0; j<nen; j++){
            connectivity->InsertCellPoint(elements->getConn(i, j));
       }     
       //printf("Rank %d:\t element %d: %d %d %d\n", myRank, i, nodes->getRemoteNode(elements->getConn(i, 0)), nodes->getRemoteNode(elements->getConn(i, 1)), nodes->getRemoteNode(elements->getConn(i, 2))); 
    }

    // Scalar property
    vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetName("Temperature");
    for(int i=0; i<nnc + nna; i++)
    {
        scalar->InsertNextValue(nodes -> getT(i));
    }
    // Previously collected data which are outputPoints, outputCells, scalarProperty, are written to
    // vtkUnstructuredGrid type grid.
    vtkSmartPointer<vtkUnstructuredGrid> unsGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unsGrid->SetPoints(outputPoints);
    unsGrid->SetCells(5,connectivity);
    unsGrid->GetPointData()->SetScalars(scalar);

    //Whatever collected in unstructured grid above is written to the "Title_myRank.vtu" file below.
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    dummy = settings->getTitle();
    dummy.append(tStep);
    dummy.append("_");
    dummy.append(strMype);
    dummy.append(".vtu");
  //  cout << dummy << endl;
    writer->SetFileName(dummy.c_str());
    writer->SetInput(unsGrid);
    writer->Write();

    // Now we write the "Title.pvtu" file which contains the informtaiton about other files.
    if(myRank == 0)
    {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = 
            vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        dummy = settings->getTitle();
        dummy.append(tStep);
        dummy.append(".pvtu");
        pwriter->SetFileName(dummy.c_str());
        pwriter->SetNumberOfPieces(numProcs);
        #if VTK_MAJOR_VERSION <= 5
            pwriter->SetInput(unsGrid);
        #else
            pwriter->SetInputData(unsGrid);
        #endif
        pwriter->Write();
    }
    ++count;
    return;
}
