#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <limits>
#include <time.h>

using namespace std;

const double PI = 3.14159265;   // Pi constant

const size_t xsd = 0;
const size_t ysd = 1;

const size_t nsd    = 2;         // number of space dimensions
const size_t nenTri = 3;         // number of element nodes
const size_t nefTri = 3;         // number of element faces

const int edgeNodesTri[3][2] = {{0,1},{1,2},{2,0}};

#endif /* CONSTANTS_H_ */
