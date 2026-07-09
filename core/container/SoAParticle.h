#pragma once

#include "ConfigVariables.hpp"
#include "math/Math.h"

struct Particles{
    std::vector<Real> mass;
    ScalarField density;

    std::vector<Real> rx;
    std::vector<Real> ry;
    std::vector<Real> rz;

    std::vector<Real> vx;
    std::vector<Real> vy;
    std::vector<Real> vz;

    std::vector<Real> fx;
    std::vector<Real> fy;
    std::vector<Real> fz;

    std::vector<uint> particle_id;
    std::vector<uint> cell_id;
};