#include "Species.h"

#include "cmath"


Species::Species(Parameters& params)
{
    dt = params.dt;
    q = params.q;
    m = params.m;
    n_particle = params.n_particle;

    charge_over_mass = q / m;

    x = new double[n_particle];
    vx = new double[n_particle];
    vy = new double[n_particle];
    vz = new double[n_particle];

    for(int i = 0; i < n_particle; i++)
    {
        double vt = sqrt(2.0 * params.temperature * const_e / m);
        double x1;
        double x2;

        x[i] = ((double)rand() / RAND_MAX) * params.lx;

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        vx[i] = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        vy[i] = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        vz[i] = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);
    }
}

void Species::dynamics(Parameters& params, double* Ex, double* rho)
{
    double Ex_local;
    double xjn, xjmxi;
    int ip;

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        /*
        xjn = x[i_particle] * params.dx_inv;
        ip  = floor(xjn);
        xjmxi = xjn - (double)ip;
        Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
        
        vx[i_particle] += charge_over_mass * Ex_local * params.dt;
        x[i_particle] += vx[i_particle] * params.dt;
        */

        interpolator1D1Order(Ex, i_particle, Ex_local);
        
        pusherBoris_test(i_particle, Ex_local);

    }

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        if(x[i_particle] < 0.0)
        {
            x[i_particle] = -x[i_particle];
        }
        else if(x[i_particle] > params.lx)
        {
            x[i_particle] = 2.0 * params.lx - x[i_particle];
        }
    }
    
    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    {
        /*
        xjn = x[i_particle] * params.dx_inv;
        ip  = floor(xjn);
        xjmxi = xjn - (double)ip;
        rho[ip] += xjmxi;
        rho[ip+1] += xjn;
        */

        projector1D1Order(rho, i_particle);
    }
}


void Species::interpolator1D1Order(double *Ex, int ipart, double &Ex_local)
{
    double dx_inv;
    double xjn, xjmxi;
    int ip;

    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
}

void Species::projector1D1Order(double *rho, int ipart)
{
    double dx_inv;
    double xjn, xjmxi;
    int ip;

    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    rho[ip] += xjmxi;
    rho[ip+1] += xjn;
} 


void Species::pusherBoris_test(int ipart, double &Ex_local)
{
    vx[ipart] += charge_over_mass * Ex_local * dt;
    x[ipart] += vx[ipart] * dt;
}
