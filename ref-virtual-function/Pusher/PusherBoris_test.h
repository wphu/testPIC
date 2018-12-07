#ifndef PUSHERBORIS_TEST_H
#define PUSHERBORIS_TEST_H

class Parameters;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_test
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_test {
public:
    //! Creator for Pusher
    PusherBoris_test(Parameters& params);
    ~PusherBoris_test();
    //! Overloading of () operator
    virtual void operator() (double *x, double *vx, int ipart, double &Ex_local);

    double dt;
    //! \todo Move mass_ in Particles_
    // mass_ relative to Species but used in the particle pusher
    double mass;
    double one_over_mass;
    double charge_over_mass;
};

#endif
