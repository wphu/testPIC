#ifndef INTERPOLATOR1D1ORDER_H
#define INTERPOLATOR1D1ORDER_H


#include "Interpolator1D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D1Order : public Interpolator1D
{

public:
    Interpolator1D1Order(Parameters& params);
    ~Interpolator1D1Order(){};

    void operator() (double *Ex, double *x, int ipart, double &Ex_local);

    inline double compute( ) {

    };

private:
    double xjn, xjmxi;
    int ip;


};//END class

#endif
