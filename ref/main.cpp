#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>

#include "Timer.h"

using namespace std;

int n_particle_per_cell = 100;
int n_space = 2000;
int nx = n_space+1;
int n_particle = n_particle_per_cell * n_space;
int n_timestep = 1000 * 2;
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

double *f;
double *a;
double *b;
double *c;
double *d;
double *e;

void solve_TDMA(double* rho, double* phi);

int main(int argc, char* argv[])
{
    double *x;
    double *vx;
    double *vy;
    double *vz;

    double *charge_density;
    double *phi;
    double *rho;
    double *Ex;

    vector<Timer> timers;
    timers.resize(5);
    timers[0].init("total time");
    timers[1].init("interpolation and push");
    timers[2].init("boundary");
    timers[3].init("projection");
    timers[4].init("solver");

    x = new double[n_particle];
    vx = new double[n_particle];
    vy = new double[n_particle];
    vz = new double[n_particle];

    charge_density = new  double[nx];
    phi = new  double[nx];
    rho = new  double[nx];
    Ex  = new  double[nx];

    f  = new  double[nx];
    a  = new  double[nx];
    b  = new  double[nx];
    c  = new  double[nx];
    d  = new  double[nx];
    e  = new  double[nx];

    for(int i = 0; i < n_particle; i++)
    {
        double vt = sqrt(2.0 * temperature * const_e / m);
        double x1;
        double x2;

        x[i] = ((double)rand() / RAND_MAX) * lx;

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

    double Ex_local;
    double xjn, xjmxi;
    int ip;

    timers[0].restart();
    for(int i_time = 0; i_time < n_timestep; i_time++)
    {
        timers[1].restart();
        for(int i_particle = 0; i_particle < n_particle; i_particle++)
        { 
            xjn = x[i_particle] * dx_inv;
            ip  = floor(xjn);
            xjmxi = xjn - (double)ip;
            Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
            
            vx[i_particle] += charge_over_mass * Ex_local * dt;
            x[i_particle] += vx[i_particle] * dt;

            /*
            if(x[i_particle] < 0.0)
            {
                x[i_particle] = -x[i_particle];
            }
            else if(x[i_particle] > lx)
            {
                x[i_particle] =2.0 * lx - x[i_particle];
            }
            */


            /*
            xjn = x[i_particle] * dx_inv;
            ip  = floor(xjn);
            xjmxi = xjn - (double)ip;
            rho[ip] += xjmxi;
            rho[ip+1] += xjn;
            */
        }
        timers[1].update();

        
        timers[2].restart();
        for(int i_particle = 0; i_particle < n_particle; i_particle++)
        { 
            if(x[i_particle] < 0.0)
            {
                x[i_particle] = -x[i_particle];
            }
            else if(x[i_particle] > lx)
            {
                x[i_particle] =2.0 * lx - x[i_particle];
            }
        }
        timers[2].update();
        

        timers[3].restart();
        for(int i_particle = 0; i_particle < n_particle; i_particle++)
        {
            xjn = x[i_particle] * dx_inv;
            ip  = floor(xjn);
            xjmxi = xjn - (double)ip;
            rho[ip] += xjmxi;
            rho[ip+1] += xjn;
        }
        timers[3].update();
        

        timers[4].restart();
        solve_TDMA(rho, phi);
        timers[4].update();
    }
    timers[0].update();


    for(int i = 0; i < timers.size(); i++)
    {
        timers[i].print();
    }
    
}

void solve_TDMA(double* rho, double* phi)
{
    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1; i++)
    {
        f[i] = -dx_square * const_ephi0_inv * rho[i];
    }

    f[0] = 0.0;
    f[nx-1] = 0.0;


    e[0] = c[0] / b[0];
    d[0] = f[0] / b[0];
    for(int i =1; i < nx-1; i++)
    {
        e[i] = c[i] / ( b[i] - a[i] * e[i-1] );
        d[i] = ( f[i] -a[i] * d[i-1] ) / ( b[i] - a[i] * e[i-1] );
    }

    phi[nx-1] = ( f[nx-1] - a[nx-1] * d[nx-2] ) / ( b[nx-1] - a[nx-1] * e[nx-2] );
    for(int i = nx-2; i >= 0; i--)
    {
        phi[i] = d[i] - e[i] * phi[i+1];
    }

}