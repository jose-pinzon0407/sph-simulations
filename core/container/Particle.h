#pragma once

#include "math/Math.h"

struct Particle{

    Particle()=default;
    Particle(Real dens, Real ms, Vec3d pos, Vec3d vel, Vec3d force);

    Real GetRho();
    Real GetMass();
    Vec3d& GetR();
    Vec3d& GetV();
    Vec3d& GetF();

    Real rho;
    Real mass=1.0;
    Vec3d r;
    Vec3d v;
    Vec3d f;

};
