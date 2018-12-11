#ifndef INTERPOLATORFACTORY_H
#define INTERPOLATORFACTORY_H

#include "Interpolator.h"
#include "Interpolator1D1Order.h"
#include "Interpolator1D2Order.h"

#include <cstdlib>

class InterpolatorFactory {
public:
    static Interpolator* create(Parameters& params) 
    {
        Interpolator* Interp = NULL;

        Interp = new Interpolator1D1Order(params);

        return Interp;
    }
};

#endif
