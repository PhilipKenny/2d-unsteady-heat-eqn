#include "setting.h"

/***************************************************************************************************
 * void setting::setting()
 * Default constructor
 **************************************************************************************************/
setting::setting()
{
    // Default values for the input parameters
    title = "Default Title";
    minfFile = "minf";
    mxyzFile = "mxyz";
    mienFile = "mien";
    mrngFile = "mrng";
//  mprmFile = "mprm";
    dataFile = "data";
    initT = 0.0;
    D = 1.0;
    source = 0.0;
    nIter = 1;
    dt = 1.0;
    restart = false;
    tolerance = 1e-6;
    numProcs = 1;
    myRank = 0;
}


/***************************************************************************************************
 * void setting::readSettings()
 * This routine calls private methods to read setting file and minf file.
 **************************************************************************************************/
void setting::readSettings()
{
    readSettingsFile();
    readMinfFile();
    
    return;
}


/***************************************************************************************************
 * void setting::readSettingsFile()
 ***************************************************************************************************
 * settings.in file contains the necessary input parameters.
 * ************************************************************************************************/
void setting::readSettingsFile()
{
    string lineString;
    string dummyString;
    char    dummyChar;
    ifstream inputFile;

   // cout << "====== Settings =============================================================" << endl;

    inputFile.open("settings.in",ios::in);
    if (inputFile.is_open()==false)
    {
        cout << "Unable to open input file! Aborting... " << endl;
        exit(0);
    }

    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "minf")
                iss >> minfFile;
            else if(dummyString == "mxyz")
                iss >> mxyzFile;
            else if(dummyString == "mien")
                iss >> mienFile;
            else if(dummyString == "mrng")
                iss >> mrngFile;
//	    else if(dummyString == "mprm")
//		iss >> mprmFile;
            else if(dummyString == "data")
                iss >> dataFile;
            else if(dummyString == "init")
                iss >> initT;
            else if(dummyString == "D")
                iss >> D;
            else if(dummyString == "S")
                iss >> source;
            else if(dummyString == "iter")
                iss >> nIter;
            else if(dummyString == "dt")
                iss >> dt;
            else if(dummyString == "restart")
                iss >> restart;
            else if(dummyString == "tolerance")
                iss >> tolerance;
            else if(dummyString == "nfg")
            {
                iss >> nFG;
                BC = new bndc [nFG+1];
            }
            else if(dummyString == "fg")
            {
                int fgnumber;
                for (int i=0; i<nFG; i++)
                {
                    iss >> fgnumber;
                    iss >> BC[fgnumber].BCType;
                    iss >> BC[fgnumber].BCValue1;
                    if (BC[fgnumber].BCType == 3)
                        iss >> BC[fgnumber].BCValue2;
                }
            }
            else
            {
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                exit(0);
            }
        }
    }

    /* cout.precision(3);
    cout << scientific;
//  cout << fixed;
    cout << "Title of the simulation                 : " << title << endl;
    cout << "Name of the minf file                   : " << minfFile << endl;
    cout << "Name of the mxyz file                   : " << mxyzFile << endl;
    cout << "Name of the mien file                   : " << mienFile << endl;
    cout << "Name of the mrng file                   : " << mrngFile << endl;
    cout << "Name of the initial distribution file   : " << dataFile << endl;
    cout << "Initial value of the dependent variable : " << initT << endl;
    cout << "Diffusion coefficient                   : " << D << endl;
    cout << "Source term                             : " << source << endl;
    cout << "Maximum number of time steps            : " << nIter << endl;
    cout << "Convergence required in max norm        : " << tolerance << endl;
    cout << "Time step size                          : " << dt << endl;
    cout << "                                         BCType\tBCValue1\tBCValue2" << endl;
    
    for (int i=1; i<nFG+1; i++)
        cout << "BC Type and values of face group " << i << "      : "
             << BC[i].BCType << '\t' << BC[i].BCValue1 << '\t' << BC[i].BCValue2 << endl;
    cout << endl;
    */

    inputFile.close();

    return;
}


/***************************************************************************************************
 * void setting::readMinf
 ***************************************************************************************************
 * minf file contains number of elements and nodes in the mesh.
 **************************************************************************************************/
void setting::readMinfFile()
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files
    double      dummyDouble;    // temperory var used for double values read from files

 //   cout << "====== Mesh =================================================================" << endl;
    file.open(minfFile.c_str(), ios::in);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    file >> dummy >> ne;
   // cout <<  "> Number of mesh elements : " << ne << endl;
    file >> dummy >> nn;
   // cout << "> Number of nodes : " << nn << endl;
   // cout << "> File read is completed: " << minfFile << endl;
    
    file.close();

    return;
}
