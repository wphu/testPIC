#include "PusherBoris_test.h"
#include "Parameters.h"

#include <iostream>
#include <cmath>

using namespace std;

PusherBoris_test::PusherBoris_test(Parameters& params)
{
    mass                = params.m;
    one_over_mass       = 1.0 / mass;
    charge_over_mass    = params.q * one_over_mass;
    dt                  = params.dt;
}

PusherBoris_test::~PusherBoris_test()
{
}

/***********************************************************************
	Only electric field -- leap-frog (Boris) scheme
***********************************************************************/
void PusherBoris_test::operator() (double *x, double *vx, int ipart, double &Ex_local)
{
    vx[ipart] += charge_over_mass * Ex_local * dt;
    x[ipart] += vx[ipart] * dt;
}
