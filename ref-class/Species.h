#ifndef SPECIES_H
#define SPECIES_H

#include "Parameters.h"

class Species
{
public:
    Species(Parameters& params);

    void dynamics(Parameters& params, double* Ex, double* rho);
    void interpolator1D1Order(double *Ex, double *x, int ipart, double &Ex_local);
    void projector1D1Order(double *rho, double *x, int ipart);
    void pusherBoris_test(double *x, double *vx, int ipart, double &Ex_local);

    double *x;
    double *vx;
    double *vy;
    double *vz;

    int n_particle;

private:
    double dt;
    
    double q;
    double m;
    double const_e = 1.6021766208e-19;
    double charge_over_mass;
};


#endif