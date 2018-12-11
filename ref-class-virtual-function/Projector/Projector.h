#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "Parameters.h"

//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector {

public:
    //! Creator for the Projector
    Projector(Parameters& params) {};
    virtual ~Projector() {};

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    virtual void operator() (double *rho, double *x, int ipart) = 0;


private:

};

#endif
