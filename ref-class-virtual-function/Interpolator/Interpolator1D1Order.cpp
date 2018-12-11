#include "Interpolator1D1Order.h"
#include "Parameters.h"

#include <cmath>
#include <iostream>

using namespace std;

Interpolator1D1Order::Interpolator1D1Order(Parameters& params) : Interpolator1D(params) {
    dx_inv = 1.0/params.dx;
}

// ---------------------------------------------------------------------------------------------------------------------
// 
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D1Order::operator() (double *Ex, double *x, int ipart, double &Ex_local)
{
    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;

}


