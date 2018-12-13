#include "Species.h"

#include "cmath"
#include <algorithm>

using namespace std;

Species::Species(Parameters& params)
{
    dt = params.dt;
    dts2 = dt / 2.0;
    q = params.q;
    m = params.m;
    n_particle = params.n_particle;

    charge_over_mass = q / m;

    timestep_sort = 100;
    i_timestep_sort = 0;

    particles.resize(n_particle);


    for(int i = 0; i < n_particle; i++)
    {
        double vt = sqrt(2.0 * params.temperature * const_e / m);
        double x1;
        double x2;

        particles[i].x = ((double)rand() / RAND_MAX) * params.lx;

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        particles[i].vx = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        particles[i].vy = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

        do 
        {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        particles[i].vz = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);
    }

    sort(particles.begin(), particles.end(), comp_particle_x);
}

void Species::dynamics(Parameters& params, double* Ex, double* rho)
{
    double Ex_local;
    double xjn, xjmxi;
    int ip;
    LocalFields E_local, B_local;

    E_local.y = 0.0;
    E_local.z = 0.0;
    B_local.x = 1.0;
    B_local.y = 1.0;
    B_local.z = 1.0;

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        interpolator1D1Order(Ex,i_particle, E_local);
        
        //pusherBoris(i_particle, E_local);
        pusherBoris(i_particle, E_local, B_local);
    }

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        if(particles[i_particle].x < 0.0)
        {
            particles[i_particle].x = -particles[i_particle].x;
        }
        else if(particles[i_particle].x > params.lx)
        {
            particles[i_particle].x = 2.0 * params.lx - particles[i_particle].x;
        }
    }

    /*
    if(i_timestep_sort == timestep_sort)
    {
        sort(particles.begin(), particles.end(), comp_particle_x);
        i_timestep_sort = 0;
    }
    else
    {
        i_timestep_sort++;
    }
    */
    
    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    {
        projector1D1Order(rho, i_particle);
    }
}


void Species::interpolator1D1Order(double *Ex, int ipart, LocalFields &E_local)
{
    double dx_inv;
    double xjn, xjmxi;
    int ip;

    xjn = particles[ipart].x * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    E_local.x = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
}

void Species::projector1D1Order(double *rho, int ipart)
{
    double dx_inv;
    double xjn, xjmxi;
    int ip;

    xjn = particles[ipart].x * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    rho[ip] += xjmxi;
    rho[ip+1] += xjn;
} 


void Species::pusherBoris(int ipart, LocalFields &E_local)
{
    particles[ipart].vx += charge_over_mass * E_local.x * dt;
    particles[ipart].x  += particles[ipart].vx * dt;
}

void Species::pusherBoris(int ipart, LocalFields &E_local, LocalFields &B_local)
{
    //double charge_over_mass_ = static_cast<double>(particles.charge(ipart))*one_over_mass_;
    double umx, umy, umz, upx, upy, upz, pxdot, pydot, pzdot;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2, Sx, Sy, Sz;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double dl;

    // --------------------------------------
    // SOLVE THE PARTICLE EQUATION OF MOTIONS
    // --------------------------------------

    // Half-acceleration in the electric field
    umx = particles[ipart].vx + charge_over_mass * E_local.x * dts2;
    umy = particles[ipart].vy + charge_over_mass * E_local.y * dts2;
    umz = particles[ipart].vz + charge_over_mass * E_local.z * dts2;

    // Rotation in the magnetic field
    alpha = charge_over_mass * dts2;
    Tx    = alpha * B_local.x;
    Ty    = alpha * B_local.y;
    Tz    = alpha * B_local.z;
    Tx2   = Tx*Tx;
    Ty2   = Ty*Ty;
    Tz2   = Tz*Tz;
    inv_det_T = 2.0/(1.0+Tx2+Ty2+Tz2);
    Sx    = Tx * inv_det_T;
    Sy    = Ty * inv_det_T;
    Sz    = Tz * inv_det_T;

    pxdot = umx + umy * Tz - umz * Ty;
    pydot = umy + umz * Tx - umx * Tz;
    pzdot = umz + umx * Ty - umy * Tx;

    upx = umx + pydot * Sz - pzdot * Sy;
    upy = umy + pzdot * Sx - pxdot * Sz;
    upz = umz + pxdot * Sy - pydot * Sx;

    // Half-acceleration in the electric field
    pxsm = upx + charge_over_mass * E_local.x * dts2;
    pysm = upy + charge_over_mass * E_local.y * dts2;
    pzsm = upz + charge_over_mass * E_local.z * dts2;

    particles[ipart].vx = pxsm;
    particles[ipart].vy = pysm;
    particles[ipart].vz = pzsm;

    // Move the particle
    particles[ipart].x += dt * particles[ipart].vx;
}
