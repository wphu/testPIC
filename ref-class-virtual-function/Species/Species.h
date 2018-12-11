#ifndef SPECIES_H
#define SPECIES_H

#include "Parameters.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "PusherFactory.h"
#include "Interpolator1D1Order_test.h"
#include "Projector1D1Order_test.h"
#include "PusherBoris_test.h"

class Species
{
public:
    Species(Parameters& params);

    void dynamics(Parameters& params, double* Ex, double* rho);

    double *x;
    double *vx;
    double *vy;
    double *vz;

    int n_particle;

private:
    double q;
    double m;
    double const_e = 1.6021766208e-19;
    double charge_over_mass;

    //Interpolator* interp;
    //Projector* proj;
    //Pusher* push;
    
    Interpolator1D1Order_test *interp;
    Projector1D1Order_test *proj;
    PusherBoris_test *push;

};


#endif