#ifndef PUSHERBORIS_H
#define PUSHERBORIS_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris : public Pusher {
public:
    //! Creator for Pusher
    PusherBoris(Parameters& params);
    ~PusherBoris();
    //! Overloading of () operator
    virtual void operator() (double *x, double *vx, int ipart, double &Ex_local);
};

#endif
