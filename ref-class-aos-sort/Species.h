#ifndef SPECIES_H
#define SPECIES_H

#include "Parameters.h"

#include <vector>

using namespace std;

struct Particle
{
    double x;
    double vx;
    double vy;
    double vz;  
};

static bool comp_particle_x (Particle a,Particle b) { return (a.x<b.x); };

// Structure containing the fields at a given position (e.g. at a Particle position)
struct LocalFields
{
    //! value of the field component along the x-direction
    double x;
    //! value of the field component along the y-direction
    double y;
    //! value of the field component along the z-direction
    double z;

};


class Species
{
public:
    Species(Parameters& params);

    void dynamics(Parameters& params, double* Ex, double* rho);
    void interpolator1D1Order(double *Ex, int ipart, LocalFields &E_local);
    void projector1D1Order(double *rho, int ipart);
    void pusherBoris(int ipart, LocalFields &E_local);
    void pusherBoris(int ipart, LocalFields &E_local, LocalFields &B_local);

    vector<Particle> particles;

    int n_particle;

private:
    double q;
    double m;
    double const_e = 1.6021766208e-19;
    double charge_over_mass;

    double dt;
    double dts2;
    int timestep_sort;
    int i_timestep_sort;
};


#endif