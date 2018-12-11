#include "Interpolator1D2Order.h"
#include "Parameters.h"

#include <cmath>
#include <iostream>

using namespace std;

Interpolator1D2Order::Interpolator1D2Order(Parameters& params) : Interpolator1D(params) {
    dx_inv = 1.0/params.dx;
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (double *Ex, double *x, int ipart, double &Ex_local)
{
    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;

}//END Interpolator1D2Order
