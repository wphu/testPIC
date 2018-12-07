#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order : public Interpolator1D
{

public:
    Interpolator1D2Order(Parameters& params);
    ~Interpolator1D2Order(){};

    void operator() (double *Ex, double *x, int ipart, double &Ex_local);


private:
    double xjn, xjmxi;
    int ip;

};//END class

#endif
