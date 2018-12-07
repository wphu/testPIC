#ifndef PROJECTOR1D_H
#define PROJECTOR1D_H

#include "Projector.h"

//----------------------------------------------------------------------------------------------------------------------
//! class Projector1D: defines a virtual method for projection in 1d3v simulations
//----------------------------------------------------------------------------------------------------------------------
class Projector1D : public Projector
{

public:
    //! Constructor for Projector1D
    Projector1D(Parameters& params) : Projector(params) {};
    virtual ~Projector1D() {};

protected:
    //! Inverse of the spatial step 1/dx
    double dx_inv;
    double dx_ov_dt;
};

#endif

