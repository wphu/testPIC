#ifndef SPECIES_H
#define SPECIES_H

#include "Parameters.h"

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
};


#endif