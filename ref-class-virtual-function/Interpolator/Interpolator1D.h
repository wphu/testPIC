#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H

#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 1D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D : public Interpolator
{
public:
    Interpolator1D(Parameters& params): Interpolator(params) {};
    virtual ~Interpolator1D() {};

    virtual void operator() (double *Ex, double *x, int ipart, double &Ex_local) = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv;
};

#endif
