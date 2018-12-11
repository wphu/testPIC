#include "PusherBoris.h"

#include <iostream>
#include <cmath>

using namespace std;

PusherBoris::PusherBoris(Parameters& params)
    : Pusher(params)
{
}

PusherBoris::~PusherBoris()
{
}

/***********************************************************************
	Only electric field -- leap-frog (Boris) scheme
***********************************************************************/
void PusherBoris::operator() (double *x, double *vx, int ipart, double &Ex_local)
{
    vx[ipart] += charge_over_mass * Ex_local * dt;
    x[ipart] += vx[ipart] * dt;
}
