#include "Interpolator1D1Order_test.h"
#include "Parameters.h"

#include <cmath>
#include <iostream>

using namespace std;

Interpolator1D1Order_test::Interpolator1D1Order_test(Parameters& params)
{
    dx_inv = 1.0/params.dx;
}

// ---------------------------------------------------------------------------------------------------------------------
// 1nd Order Interpolation of the fields at a the particle position (2 nodes are used)
// only has Electric field, no Magnetic field
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D1Order_test::operator() (double *Ex, double *x, int ipart, double &Ex_local)
{
    xjn = x[ipart] * dx_inv;
    ip  = floor(xjn);
    xjmxi = xjn - (double)ip;
    Ex_local = Ex[ip] * xjmxi + Ex[ip+1] * xjn;
}



