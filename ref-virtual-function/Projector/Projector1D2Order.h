#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
    Projector1D2Order(Parameters& params);
    ~Projector1D2Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    void operator() (double *rho, double *x, int ipart);

private:
    double dx_ov_dt;
};

#endif
