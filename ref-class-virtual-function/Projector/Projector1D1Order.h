#ifndef PROJECTOR1D1ORDER_H
#define PROJECTOR1D1ORDER_H

#include "Projector1D.h"

class Projector1D1Order : public Projector1D {
public:
    Projector1D1Order(Parameters& params);
    ~Projector1D1Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    void operator() (double *rho, double *x, int ipart);


private:
    double xjn, xjmxi;
    int ip;
};

#endif
