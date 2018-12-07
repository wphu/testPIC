#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

//#include "Parameters.h"

class Parameters;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator
{
public:
    Interpolator(Parameters& params) {};
    virtual ~Interpolator() {};

    virtual void operator() (double *Ex, double *x, int ipart, double &Ex_local) = 0;

private:

};//END class

#endif
