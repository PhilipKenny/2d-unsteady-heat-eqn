#ifndef POSTPROCESSOR_H_
#define POSTPROCESSOR_H_

#include "solver.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>

class postProcessor
{
    private:
        //VARIABLES
        setting*  settings;   // a local pointer for the settings
        element* elements;
        node* nodes;
        double          minT;       // min value of the Temperature field
        double          maxT;
        
        //METHODS
        void evaluateLimits();
        void vtkVisualization();
    protected:

    public:
        void postProcessorControl(setting* argSettings, node* argNodes, element* argElements);
};


#endif /* POSTPROCESSOR_H_ */
