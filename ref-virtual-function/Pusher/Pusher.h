#ifndef PUSHER_H
#define PUSHER_H

#include "Parameters.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class Pusher
{

public:
    //! Creator for Pusher
    Pusher(Parameters& params);
    virtual ~Pusher();

    //! Overloading of () operator
    virtual void operator() (double *x, double *vx, int ipart, double &Ex_local) = 0;

protected:
    double dt;
    //! \todo Move mass_ in Particles_
    // mass_ relative to Species but used in the particle pusher
    double mass;
    double one_over_mass;
    double charge_over_mass;



};//END class

#endif
