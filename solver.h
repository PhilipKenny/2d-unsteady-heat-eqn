#ifndef SOLVER_H_
#define SOLVER_H_

#include "node.h"
#include "element.h"
#include "postProcessor.h"

/***************************************************************************************************
 * Solver class
***************************************************************************************************/
class femSolver
{
    private:
        setting*    stngs;  // a local pointer for the settings
        node*       nodes;  // a local pointer for the nodes
        element*    elems;  // a local pointer for the elements

        // PRIVATE METHODS
        bool explicitSolver(const double& tolerance);

    protected:

    public:
        // CONSTRUCTOR //
        femSolver(setting* argSet, node* argNode, element* argElem)
        {
            stngs = argSet;
            nodes = argNode;
            elems = argElem;
        };

        // DESTRUCTOR
        ~femSolver(){};
        
        // INTERFACE FUNCTION
        void solve();

};

#endif /* SOLVER_H_ */


