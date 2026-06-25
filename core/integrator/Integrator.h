#pragma once

/***
 * A simple euler integrator as a default
 */

#include "math/Math.h"
#include "Domain.h"
#include "container/Container.h"


class Integrator{
    public:

    Integrator(Real delta_t, Domain* domain);

    void UpdatePositions(Container& particles);
    void UpdateVelocities(Container& particles);


    protected:
    bool ApplyBoundaryConditions(Vec3d& r, Vec3d& v);
    // bool ApplyBoundaryConditions(Particle& particle);

    protected:
    Domain* domain=nullptr;
    double dt=0.0;
};