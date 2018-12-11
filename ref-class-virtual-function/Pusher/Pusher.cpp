#include "Pusher.h"

Pusher::Pusher(Parameters& params)
{
    mass                = params.m;
    one_over_mass       = 1.0 / mass;
    charge_over_mass    = params.q * one_over_mass;
    dt                  = params.dt;
}

Pusher::~Pusher()
{
}
