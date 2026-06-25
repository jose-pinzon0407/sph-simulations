#pragma once

#include "math/Math.h"

class Container;
struct Particle{
    friend class Container;
    public:
    Particle()=default;
    Particle(Real dens, Real ms, Vec3d& pos, Vec3d& vel, Vec3d& force);

    Real GetRho();
    Real GetMass();
    Vec3d& GetR();
    Vec3d& GetV();
    Vec3d& GetF();

    public:
    Real rho;
    Real mass=1.0;
    Vec3d r;
    Vec3d v;
    Vec3d f;

};
