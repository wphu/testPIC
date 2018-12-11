#include "Species.h"

#include "cmath"


Species::Species(Parameters& params)
{
    dt = params.dt;
    dts2 = dt / 2.0;
    q = params.q;
    m = params.m;
    n_particle = params.n_particle;

    charge_over_mass = q / m;

    x.resize(n_particle);
    vx.resize(n_particle);
    vy.resize(n_particle);
    vz.resize(n_particle);

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
    LocalFields Epart_local, Bpart_local;

    Epart_local.y = 0.0;
    Epart_local.z = 0.0;
    Bpart_local.x = 1.0;
    Bpart_local.y = 1.0;
    Bpart_local.z = 1.0;

    for(int i_particle = 0; i_particle < n_particle; i_particle++)
    { 
        interpolator1D1Order(Ex, i_particle, Epart_local);
        
        //pusherBoris(i_particle, Epart_local);
        pusherBoris(i_particle, Epart_local, Bpart_local);

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
        projector1D1Order(rho, i_particle);
    }
}


void Species::interpolator1D1Order(double *Ex, int ipart, LocalFields &Epart)
{
    double dx_inv;
    double xjn, xjmxi;
    int ip;

    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    Epart.x = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
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


void Species::pusherBoris(int ipart, LocalFields &Epart)
{
    vx[ipart] += charge_over_mass * Epart.x * dt;
    x[ipart] += vx[ipart] * dt;
}

void Species::pusherBoris(int ipart, LocalFields &Epart, LocalFields &Bpart)
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
    umx = vx[ipart] + charge_over_mass * Epart.x * dts2;
    umy = vy[ipart] + charge_over_mass * Epart.y * dts2;
    umz = vz[ipart] + charge_over_mass * Epart.z * dts2;

    // Rotation in the magnetic field
    alpha = charge_over_mass * dts2;
    Tx    = alpha * Bpart.x;
    Ty    = alpha * Bpart.y;
    Tz    = alpha * Bpart.z;
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
    pxsm = upx + charge_over_mass * Epart.x * dts2;
    pysm = upy + charge_over_mass * Epart.y * dts2;
    pzsm = upz + charge_over_mass * Epart.z * dts2;

    vx[ipart] = pxsm;
    vy[ipart] = pysm;
    vz[ipart] = pzsm;

    // Move the particle
    x[ipart] += dt * vx[ipart];

}
