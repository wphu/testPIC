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
    void interpolator1D1Order(double *Ex, int i_space, int i_particle, LocalFields &E_local);
    void projector1D1Order(double *rho, int i_space, int i_particle);
    void pusherBoris(int i_space, int i_particle, LocalFields &E_local);
    void pusherBoris(int i_space, int i_particle, LocalFields &E_local, LocalFields &B_local);
    void sort_particles();
    void sort_particles_new();

    vector< vector<Particle> > particles;
    vector< vector<int> > list_removed_particles;

    int n_particle;

private:
    double q;
    double m;
    double const_e = 1.6021766208e-19;
    double charge_over_mass;

    double dt;
    double dts2;
    double dx;
    double dx_inv;
    double lx;
    int timestep_sort;
    int i_timestep_sort;
};


#endif