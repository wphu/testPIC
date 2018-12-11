#include "Projector1D1Order.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector1D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector1D1Order::Projector1D1Order (Parameters& params) : Projector1D(params)
{
    dx_inv  = 1.0/params.dx;
    dx_ov_dt = params.dx / params.dt;
}

Projector1D1Order::~Projector1D1Order()
{
}



// ---------------------------------------------------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D1Order::operator() (double *rho, double *x, int ipart)
{
    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    rho[ip] += xjmxi;
    rho[ip+1] += xjn;
} 
