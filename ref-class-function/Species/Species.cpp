#include "Species.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "PusherFactory.h"
#include "Interpolator1D1Order_test.h"
#include "Projector1D1Order_test.h"
#include "PusherBoris_test.h"

#include "cmath"


Species::Species(Parameters& params)
{
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

    //Interpolator* interp    = InterpolatorFactory::create(params);
    //Projector* proj           = ProjectorFactory::create(params);
    //Pusher* push              = PusherFactory::create(params);
    
    Interpolator1D1Order_test *interp    = new Interpolator1D1Order_test(params);
    Projector1D1Order_test *proj         = new Projector1D1Order_test(params);
    PusherBoris_test *push               = new PusherBoris_test(params);

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        (*interp) (Ex, x, i_particle, Ex_local);
        (*push) (x, vx, i_particle, Ex_local);

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
        (*proj) (rho, x, i_particle);
    }
}