#include "gaussq.h"

/***************************************************************************************************
 * void GaussQ::setupPointsAndWeights()
 ***************************************************************************************************
 * GAUSS QUADRATURE PIONTS AND WEIGHTS
***************************************************************************************************/
void GaussQ::setupPointsAndWeights()
{
    if (nGQP == 7)
    {
        ksi[0]    = 0.333333333333333;
        eta[0]    = 0.333333333333333;
        weight[0] = 0.225 / 2.0;

        ksi[1]    = 0.059715871789770;   
        eta[1]    = 0.470142064105115;
        weight[1] = 0.132394152788 / 2.0;

        ksi[2]    = 0.470142064105115;   
        eta[2]    = 0.059715871789770;
        weight[2] = 0.132394152788 / 2.0;

        ksi[3]    = 0.470142064105115;   
        eta[3]    = 0.470142064105115;
        weight[3] = 0.132394152788 / 2.0;

        ksi[4]    = 0.101286507323456;   
        eta[4]    = 0.797426985353087;
        weight[4] = 0.125939180544 / 2.0;

        ksi[5]    = 0.101286507323456;   
        eta[5]    = 0.101286507323456;
        weight[5] = 0.125939180544 / 2.0;

        ksi[6]    = 0.797426985353087;   
        eta[6]    = 0.101286507323456;
        weight[6] = 0.125939180544 / 2.0;
    }
    else
    {
        cout << "Quadrature rule not defined." << endl;
        exit(0);
    }
    return;
}

/***************************************************************************************************
 * void GaussQ::evaluateShapeFunctions()
 ***************************************************************************************************
 * EVALUATES SHAPE FUNCTIONS FOR LINEAR TRIANGULAR ELEMENT
***************************************************************************************************/
void GaussQ::evaluateShapeFunctions()
{
    for(int iGQP=0; iGQP<nGQP; iGQP++)
    {
        if (nen == 3)
        {
            setS(iGQP, 0, 1.0-ksi[iGQP]-eta[iGQP]);
            setS(iGQP, 1, ksi[iGQP]);
            setS(iGQP, 2, eta[iGQP]);
            
            setDSdKsi(iGQP, 0, -1.0);
            setDSdKsi(iGQP, 1,  1.0);
            setDSdKsi(iGQP, 2,  0.0);

            setDSdEta(iGQP, 0, -1.0);
            setDSdEta(iGQP, 1,  0.0);
            setDSdEta(iGQP, 2,  1.0);
        }
        else
        {
            cout << "Element type not defined" << endl;
            exit(0);
        }
    }

    return;
}



