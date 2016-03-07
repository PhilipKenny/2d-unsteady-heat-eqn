#include <string>
#include "utility.h"
#include <iostream>
#include <ostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{


  if (argc < 2) return 1;

  std::string numProcs(argv[1]);
  std::string pwd(argv[2]);
  size_t ne, nn;
  int nen = 3;
  int nsd = 2;
  std::string mien_neu = string("mien.");
  std::string mxyz_neu("mxyz.");
  std::string mrng_neu("mrng.");
  std::string nprm("nprm.");
  std::string mprm("mprm.");     

 // pwd = "../mesh-Rectangle/finemesh/";

 
    fstream    file;           // file name obj
    std::string dummy = string(pwd) + "minf";          // dummy string to hold names
    file.open(dummy.c_str(), ios::in);
 //   file.open("../mesh-Rectangle/finemesh/minf", ios::in);

    file >> dummy >> ne;
    file >> dummy >> nn;
//  cout << nn << endl;
    file.close();

   double* coord = new double[nn*nsd];
   int* conn = new int[ne*nen];
   int* facer = new int[ne*nen];
   
    dummy = pwd + "mxyz";          // dummy string to hold names
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
 

    file.seekg (0, ios::beg);
    file.read((char*)coord, nsd*sizeof(double)*nn);
    utility::swapBytes((char*)coord, nsd*nn, sizeof(double));
    cout << "> File read is completed: " << dummy << endl;
    file.close();


    cout << dummy << endl; 

    dummy = pwd + "mien";          // dummy string to hold names
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);

    file.seekg (0, ios::beg);
    file.read((char*)conn, nen*sizeof(int)*ne);
    utility::swapBytes((char*)conn, nen*ne, sizeof(int));
    file.close();


    cout << dummy << endl; 

    dummy = pwd + "mrng";          // dummy string to hold names
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    cout << dummy << endl; 

    file.seekg (0, ios::beg);
    file.read((char*)facer, nen*sizeof(int)*ne);
    utility::swapBytes((char*)facer, nen*ne, sizeof(int));
    file.close(); 


    int* nperm = new int[nn]();
    int* mperm = new int[ne];

    ostringstream int2str; 
    int2str << setfill('0') << setw(5) << numProcs;
    nprm.append(int2str.str()); 
    mprm.append(int2str.str()); 
    
    dummy = pwd + nprm;
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    file.seekg (0, ios::beg);
    file.read((char*)nperm, sizeof(int)*nn);
    utility::swapBytes((char*)nperm, nn, sizeof(int));
    file.close(); 

    dummy = pwd + mprm;
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    file.seekg (0, ios::beg);
    file.read((char*)mperm, sizeof(int)*ne);
    utility::swapBytes((char*)mperm, ne, sizeof(int));
    file.close();

    for (size_t in = 0; in < nn; ++in)
        --nperm[in];
    for (size_t ie = 0; ie < ne; ++ie)
        --mperm[ie];

    

  // Swap bytes here
  // output to the file you want 
  // I propose EDOCOTWOHAEDIONEVAHUOYPILIHP.txt
    double* coordneu = new double[nn*nsd];
    int* connneu = new int[ne*nen];
    int* facerneu = new int[ne*nen];
  
//    cout << nprm << endl; 
 //   cout << sizeof(int)*nn <<endl;
    for (size_t in = 0; in < nn; ++in)
    {
        coordneu[nperm[in]*nsd] = coord[in*nsd]; 
        coordneu[nperm[in]*nsd+1] = coord[in*nsd+1]; 
//        cout << nperm[in] << endl;
    }

    for (size_t ie = 0; ie < ne; ++ie)
        for (int ien = 0; ien < nen; ++ien)
        {
            size_t val = conn[ie*nen+ien] - 1;
            connneu[ie*nen+ien] = nperm[val] + 1;
 //           facerneu[ie*nen+ien] = nperm[val] + 1;
        }

    for (size_t ie = 0; ie < ne; ++ie)
    {
        conn[mperm[ie]*nen] = connneu[ie*nen]; 
        conn[mperm[ie]*nen+1] = connneu[ie*nen+1]; 
        conn[mperm[ie]*nen+2] = connneu[ie*nen+2]; 

        facerneu[mperm[ie]*nen] = facer[ie*nen]; 
        facerneu[mperm[ie]*nen+1] = facer[ie*nen+1]; 
        facerneu[mperm[ie]*nen+2] = facer[ie*nen+2]; 
    }
  
    utility::swapBytes((char*)coordneu, nsd*nn, sizeof(double));
    dummy = pwd + mxyz_neu + int2str.str();
    ofstream fout(dummy.c_str(), ios::out|ios::binary);
    fout.write((char*)coordneu, sizeof(double)*nsd*nn);
    fout.close();

    utility::swapBytes((char*)conn, nen*ne, sizeof(int));
    dummy = pwd + mien_neu + int2str.str();
    ofstream fout2(dummy.c_str(), ios::out|ios::binary);
    fout2.write((char*)conn, sizeof(int)*ne*nen);
    fout2.close();

    utility::swapBytes((char*)facerneu, nen*ne, sizeof(int));
    dummy = pwd + mrng_neu + int2str.str();
    ofstream fout3(dummy.c_str(), ios::out|ios::binary);
    fout3.write((char*)facerneu, sizeof(int)*ne*nen);
    fout3.close();
    return 0;
}
