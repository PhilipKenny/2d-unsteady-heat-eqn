#ifndef GAUSSQ_H_
#define GAUSSQ_H_

#include "setting.h"

/***************************************************************************************************
 * QUADRATURE LEVEL DATA STRUCTURE
 ***************************************************************************************************
 * At each quarature point we define:
 *      coordinates on Master element, namely ksi and eta
 *      quadrature weights
 *      shape functions
 *      shape function derivatives wrt. ksi and eta
 *      shape function derivatives wrt. x and y
 *      Jacobian determinant
***************************************************************************************************/
class GaussQ
{
    private:
        // PRIVATE VARIABLES
        size_t ne;
        size_t nen;
        size_t nGQP;

        // PRIVATE METHODS
        void setupPointsAndWeights();
        void evaluateShapeFunctions();

    protected:
        // PROTECTED VARIABLES
        double* ksi;    // ksi for each GQ point              [nGQP]
        double* eta;    // eta for each GQ point              [nGQP]
        double* weight; // weight of each GQ point            [nGQP]
        double* S;      // Shape functions                    [nGQP x nen]
        double* dSdKsi; // ksi derivatives of shape functions [nGQP x nen]
        double* dSdEta; // eta derivatives of shape functions [nGQP x nen]
        double* detJ;   // Jacobian determinant. Dimension is [ne x nGQP]
        double* dSdx;   // Shape function derivatives wrt. x. [ne x nGQP x nen]
        double* dSdy;   // Shape function derivatives wrt. y. [ne x nGQP x nen]

        // PROTECTED METHODS
        inline void setS      (           size_t iGQP, size_t ien, double value) {S[iGQP*nen+ien] = value;};
        inline void setDSdKsi (           size_t iGQP, size_t ien, double value) {dSdKsi[iGQP*nen+ien] = value;};
        inline void setDSdEta (           size_t iGQP, size_t ien, double value) {dSdEta[iGQP*nen+ien] = value;};
        inline void setDetJ   (size_t ie, size_t iGQP,             double value) {detJ[ie*nGQP+iGQP] = value;};
        inline void setDSdx   (size_t ie, size_t iGQP, size_t ien, double value) {dSdx[ie*nGQP*nen+iGQP*nen+ien] = value;};
        inline void setDSdy   (size_t ie, size_t iGQP, size_t ien, double value) {dSdy[ie*nGQP*nen+iGQP*nen+ien] = value;};

        inline double getWeight (size_t iGQP            )  {return weight[iGQP];};
        inline double getS      (size_t iGQP, size_t ien)  {return S[iGQP*nen+ien];};
        inline double getDSdKsi (size_t iGQP, size_t ien)  {return dSdKsi[iGQP*nen+ien];};
        inline double getDSdEta (size_t iGQP, size_t ien)  {return dSdEta[iGQP*nen+ien];};
        inline double getDetJ   (size_t ie,   size_t iGQP) {return detJ[ie*nGQP+iGQP];};
        inline double getDSdx   (size_t ie,   size_t iGQP, size_t ien) {return dSdx[ie*nGQP*nen+iGQP*nen+ien];};
        inline double getDSdy   (size_t ie,   size_t iGQP, size_t ien) {return dSdy[ie*nGQP*nen+iGQP*nen+ien];};

    public:
        // DEFAULT CONSTRUCTOR //
        GaussQ(){};

        // CONSTRUCTOR //
        GaussQ(size_t argNe, size_t argNen, size_t argNGQP)
        {
            ne     = argNe;
            nen    = argNen;
            nGQP   = argNGQP;
            ksi    = new double[nGQP];
            eta    = new double[nGQP];
            weight = new double[nGQP];
            S      = new double[nGQP * nen];
            dSdKsi = new double[nGQP * nen];
            dSdEta = new double[nGQP * nen];
            detJ   = new double[ne * nGQP];
            dSdx   = new double[ne * nGQP * nen];
            dSdy   = new double[ne * nGQP * nen];

            setupPointsAndWeights();
            evaluateShapeFunctions();
        };

        // DESTRUCTOR //
        ~GaussQ()
        {
            delete[] ksi;
            delete[] eta;
            delete[] weight;
            delete[] S;
            delete[] dSdKsi;
            delete[] dSdEta;
            delete[] detJ;
            delete[] dSdx;
            delete[] dSdy;
        };
};

#endif /* GAUSSQ_H_ */
