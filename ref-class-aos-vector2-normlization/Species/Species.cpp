#include "Species.h"

#include "cmath"
#include <algorithm>
#include <iostream>

using namespace std;

Species::Species(Parameters& params)
{
    dt = params.dt;
    dts2 = dt / 2.0;
    dx = params.dx;
    dx_inv = 1.0 / dx;
    lx = params.lx;
    dt_over_dx = dt / dx;
    q = params.q;
    m = params.m;
    n_particle = params.n_particle;
    n_space = params.n_space;

    charge_over_mass = q / m;

    timestep_sort = 100;
    i_timestep_sort = 0;

    particles.resize(params.n_space);
    list_removed_particles.resize(params.n_space);
    for(int i = 0; i < params.n_space; i++)
    {
        particles[i].reserve(1.5 * params.n_particle_per_cell);
        particles[i].resize(params.n_particle_per_cell);
    }

    for(int i_space = 0; i_space < params.n_space; i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        {
            double vt = sqrt(2.0 * params.temperature * const_e / m);
            double x1;
            double x2;

            particles[i_space][i_particle].x = ((double)rand() / RAND_MAX);

            do 
            {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_space][i_particle].vx = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

            do 
            {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_space][i_particle].vy = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

            do 
            {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_space][i_particle].vz = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);
        }
    }


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

    //cout<<"1111 "<<endl;
    for(int i_space = 0; i_space < params.n_space; i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        { 
            interpolator1D1Order(Ex,i_space, i_particle, E_local);
            
            E_local.x = 0.0;
            //pusherBoris(i_space, i_particle, E_local);
            pusherBoris(i_space, i_particle, E_local, B_local);
        }
    }

    //cout<<"222 "<<endl;
    //sort_particles();
    sort_particles();

    //cout<<"3333 "<<endl;
    for(int i_space = 0; i_space < params.n_space; i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        { 
            projector1D1Order(rho, i_space, i_particle);
        }
    }
}


void Species::interpolator1D1Order(double *Ex, int i_space, int i_particle, LocalFields &E_local)
{
    double xjn;

    xjn = particles[i_space][i_particle].x;
    E_local.x = Ex[i_space] * (1.0 - xjn) + Ex[i_space+1] * xjn;
}

void Species::projector1D1Order(double *rho, int i_space, int i_particle)
{
    double xjn;

    xjn = particles[i_space][i_particle].x;
    rho[i_space] += (1.0 - xjn);
    rho[i_space+1] += xjn;
} 


void Species::pusherBoris(int i_space, int i_particle, LocalFields &E_local)
{
    particles[i_space][i_particle].vx += charge_over_mass * E_local.x * dt;
    particles[i_space][i_particle].x  += particles[i_space][i_particle].vx * dt_over_dx;
}

void Species::pusherBoris(int i_space, int i_particle, LocalFields &E_local, LocalFields &B_local)
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
    umx = particles[i_space][i_particle].vx + charge_over_mass * E_local.x * dts2;
    umy = particles[i_space][i_particle].vy + charge_over_mass * E_local.y * dts2;
    umz = particles[i_space][i_particle].vz + charge_over_mass * E_local.z * dts2;

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

    particles[i_space][i_particle].vx = pxsm;
    particles[i_space][i_particle].vy = pysm;
    particles[i_space][i_particle].vz = pzsm;

    // Move the particle
    particles[i_space][i_particle].x += dt_over_dx * particles[i_space][i_particle].vx;
}




void Species::sort_particles()
{
    int i_space_change;
    int n_space_check_boundary;

    n_space_check_boundary = 5;


    //cout<<"aaa"<<endl;
    for(int i_space = 0; i_space < n_space_check_boundary; i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        {
            i_space_change =  floor(particles[i_space][i_particle].x);

            if(i_space + i_space_change < 0)
            {
                particles[i_space][i_particle].x = -particles[i_space][i_particle].x - 2.0 * i_space;
            }
        }
    }
    //cout<<"bbb"<<endl;
    for(int i_space = particles.size() - n_space_check_boundary; i_space < particles.size(); i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        {
            i_space_change =  floor(particles[i_space][i_particle].x);

            if(i_space + i_space_change >= particles.size())
            {
                particles[i_space][i_particle].x = -particles[i_space][i_particle].x - 2.0 * (i_space - n_space);
            }
        }
    }

    //cout<<"ccc"<<endl;
    for(int i_space = 0; i_space < particles.size(); i_space++)
    {
        for(int i_particle = 0; i_particle < particles[i_space].size(); i_particle++)
        {
            i_space_change =  floor(particles[i_space][i_particle].x);

            if(i_space_change != 0)
            {
                particles[i_space][i_particle].x = abs(particles[i_space][i_particle].x - i_space_change);
                particles[i_space + i_space_change].push_back(particles[i_space][i_particle]);
                list_removed_particles[i_space].push_back(i_particle);
            }
        }
    }
    //cout<<"ddd"<<endl;
    for(int i_space = 0; i_space < particles.size(); i_space++)
    {
        int i_particle_last = particles[i_space].size() - 1;
        for(int i_list = 0; i_list < list_removed_particles[i_space].size(); i_list++)
        {
            int i_particle = list_removed_particles[i_space][list_removed_particles[i_space].size() - 1 - i_list];
            if(i_particle == i_particle_last)
            {
                i_particle_last--;
            }
            else
            {
                swap(particles[i_space][i_particle], particles[i_space][i_particle_last]);
                i_particle_last--;
            }
        }
        particles[i_space].resize(particles[i_space].size() - list_removed_particles[i_space].size());
        list_removed_particles[i_space].clear();
    } 
    //cout<<"eee"<<endl;
}