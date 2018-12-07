#include "Projector1D2Order.h"
#include "Parameters.h"

#include <cmath>
#include <iostream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector1D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector1D2Order::Projector1D2Order(Parameters& params)
    :Projector1D(params)
{
    dx_inv  = 1.0/params.dx;
    dx_ov_dt = params.dx / params.dt;

}

Projector1D2Order::~Projector1D2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2Order::operator() (double *rho, double *x, int ipart)
{



} 


