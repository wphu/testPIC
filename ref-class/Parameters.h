#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters
{
public:

Parameters()
{
    n_particle_per_cell = 100;
    n_space = 2000;
    nx = n_space+1;
    n_particle = n_particle_per_cell * n_space;
    n_timestep = 1000;
    dt = 1.0e-12;
    lx = 1.0e-2;
    temperature = 20;
    density = 1.0e19;
    m = 9.109382616e-31;
    q = -1.6021766208e-19;
    weight = density / n_particle_per_cell;

    dx = lx / n_space;
    dx_inv = 1.0 / dx;
    dx_square = dx * dx;
    charge_over_mass = q / m;

    const_ephi0 = 8.854187817620389e-12;
    const_ephi0_inv = 1.0 / const_ephi0;
    const_e = 1.6021766208e-19;
    /*
    const_e = 1.6021766208e-19;
    const_emass = 9.10938356e-31;
    const_c = 299792458.0;
    const_ephi0 = 8.854187817620389e-12;
    const_pi = 3.1415926;
    const_boltz = 1.3806e-23;
    const_h = 6.62606957e-34;
    */
}

    int n_particle_per_cell = 100;
    int n_space = 2000;
    int nx = n_space+1;
    int n_particle = n_particle_per_cell * n_space;
    int n_timestep = 1000;
    double dt = 1.0e-12;
    double lx = 1.0e-2;
    double temperature = 20;
    double density = 1.0e19;
    double m = 9.109382616e-31;
    double q = -1.6021766208e-19;
    double weight = density / n_particle_per_cell;

    double dx = lx / n_space;
    double dx_inv = 1.0 / dx;
    double dx_square = dx * dx;
    double charge_over_mass = q / m;

    double const_ephi0 = 8.854187817620389e-12;
    double const_ephi0_inv = 1.0 / const_ephi0;
    double const_e = 1.6021766208e-19;
    /*
    const_e = 1.6021766208e-19;
    const_emass = 9.10938356e-31;
    const_c = 299792458.0;
    const_ephi0 = 8.854187817620389e-12;
    const_pi = 3.1415926;
    const_boltz = 1.3806e-23;
    const_h = 6.62606957e-34;
    */
};


#endif