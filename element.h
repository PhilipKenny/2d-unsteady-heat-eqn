#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "node.h"
#include "gaussq.h"
#include <vector>

#define swap(a,b) {double temp; temp=a; a=b; b=temp;}

/***************************************************************************************************
 * ELEMENT LEVEL DATA STRUCTURE. 
 * Derived from GaussQ class
 ***************************************************************************************************
 * Each element has
 *      connectivity [ne x nen]
 *      face groups defined for each face [ne x nef]
 *      Mass matrix [ne x nen]
 *      Forcing matrix [ne x nen]
 *      Stiffness matrix [ne x nen x nen]
***************************************************************************************************/
class element : public GaussQ
{
    private:
        // PRIVATE VARIABLES
        size_t  ne;     // Number of elements
        size_t  nen;    // Number of element nodes
        size_t  nef;    // Number of element faces
        size_t  nGQP;   // Number of Gauss quadrature points
        int*    conn;   // Connectivity             [ne x nen]
        int*    FG;     // Face groups              [ne x nef]
        double* M;      // Element mass matrix      [ne * nen]
        double* F;      // Element RHS              [ne * nen]
        double* K;      // Element stiffness matrix [ne * nen * nen]
        // @Klim        // The *.c was also updated
        double* RHS;    // Nodal global RHS         [ne * nen] 
        double* bound;   // Boundary condition vals  [ne * nen]
        size_t  maxElems;         


        // PRIVATE METHODS
        void setM (size_t ie, size_t ien, double value) {M[ie*nen+ien] = value;};
        void setF (size_t ie, size_t ien, double value) {F[ie*nen+ien] = value;};
        void addF (size_t ie, size_t ien, double value) {F[ie*nen+ien] += value;};
        void setK (size_t ie, size_t ien, size_t jen, double value) {K[ie*nen*nen+ien*nen+jen] = value;};
        // @Klim
        void setRHS (size_t ie, size_t ien, double value) {RHS[ie*nen+ien] = value;};
        void addRHS (size_t ie, size_t ien, double value) {RHS[ie*nen+ien] += value;};
        void setBound (size_t ie, size_t ien, double value) {bound[ie*nen+ien] = value;};
        void addBound (size_t ie, size_t ien, double value) {bound[ie*nen+ien] += value;};
        void addK (size_t ie, size_t ien, size_t jen, double value) {K[ie*nen*nen+ien*nen+jen] += value;};

        void readMien(setting*);
        void readMrng(setting*);
        void readMienPar(setting* settings);
	void readMrngPar(setting* settings);

	void calculateJacobian(node*);
        void calculateElementMatrices(setting*, node*);
        void applyBoundaryConditions(setting*, node*);
        // My functions
        void setStiffnessElement(size_t ie, size_t ien, size_t jen, double diffuse); 
        void setForcingElement(size_t ie, size_t ien, double f);
        void setConsistentMassElement(size_t ie, size_t ien, size_t jen, double* consistentMassMatrix);
        void setLumpedMassMatrix(size_t ie, double* consistentMassMatrix);
        void printElementStiffnessMatrix(size_t ie);
        void printElementMassMatrix(size_t ie);
        void modifyConn(node* nodes);       


    protected:

    public:
        // CONSTRUCTOR //
        element(size_t argNe, size_t argNen, size_t argNef, size_t argNGQP)
        : GaussQ(argNe, argNen, argNGQP)
        {
            ne   = argNe;
            nen  = argNen;
            nef  = argNef;
            nGQP = argNGQP;
            conn = new int[ne*nen];
            FG   = new int[ne*nef];
            M    = new double[ne*nen];
            F    = new double[ne*nen];
            K    = new double[ne*nen*nen];
            // @Klim
            RHS  = new double[ne*nen];
            bound = new double[ne*nen];
        };

        // DESTRUCTOR //
        ~element()
        {
            delete[] conn;
            delete[] FG;
            delete[] K;
            delete[] F;
            delete[] M;
            delete[] RHS;
            delete[] bound;
        };

        // PUBLIC GETTERS
        size_t  getNen  ()          {return nen;};
        size_t  getNe   ()	    {return ne;};		
	size_t  getNGQP ()          {return nGQP;};
        int     getConn (size_t ie, size_t ien) {return conn[ie*nen+ien];};
        double* getMptr (size_t ie) {return &M[ie*nen];};
        double* getFptr (size_t ie) {return &F[ie*nen];};
        double* getKptr (size_t ie) {return &K[ie*nen*nen];};
        // @ Klim
        double*   getRHSptr(size_t ie) {return &RHS[ie*nen];};
        double    getRHS  (size_t ie, size_t ien) {return RHS[ie*nen+ien];};
        double    getF (size_t ie, size_t ien) {return F[ie*nen+ien];};
        int       getFG   (size_t ie, size_t ief) {return FG[ie*nef+ief];};
	void 	  calculateRHS(node* nodes, setting* settings);       
        double*   getBoundptr(size_t ie) {return &bound[ie*nen];};
        double    getBound(size_t ie, size_t ien) {return bound[ie*nen+ien];};
        size_t    getMaxElems()      {return maxElems;};

        // PUBLIC INTERFACE
        void createElements(setting*, node*);
        void calculateOffsets(node*);
};

#endif /* ELEMENT_H_ */
