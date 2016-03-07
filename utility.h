#ifndef UTILITY_H_
#define UTILITY_H_

#include "constants.h"
#include <vector>

namespace utility
{
    extern void swapBytes (char *array, int nelem, int elsize);
    extern void resizeVec (std::vector<std::vector<int> > &input, const std::vector<std::vector<int> > & pattern);
    extern void resizeVec (std::vector<std::vector<double> > &input, const std::vector<std::vector<int> > & pattern);
}

#endif /* UTILITY_H_ */
