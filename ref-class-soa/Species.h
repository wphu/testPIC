#ifndef SPECIES_H
#define SPECIES_H

#include "Parameters.h"

#include <vector>

using namespace std;


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
    void interpolator1D1Order(double *Ex, int ipart, LocalFields &Epart);
    void projector1D1Order(double *rho, int ipart);
    void pusherBoris(int ipart, LocalFields &Epart);
    void pusherBoris(int ipart, LocalFields &Epart, LocalFields &Bpart);

    vector<double> x;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;

    int n_particle;

private:
    double dt;
    double dts2;
    
    double q;
    double m;
    double const_e = 1.6021766208e-19;
    double charge_over_mass;
};


#endif