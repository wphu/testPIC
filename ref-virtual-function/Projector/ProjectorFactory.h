#ifndef PROJECTORFACTORY_H
#define PROJECTORFACTORY_H

#include "Projector.h"
#include "Projector1D1Order.h"
#include "Projector1D2Order.h"

class ProjectorFactory {
public:
    static Projector* create(Parameters& params) 
    {
        Projector* Proj = NULL;

        Proj = new Projector1D1Order(params);

        return Proj;
    }

};

#endif
