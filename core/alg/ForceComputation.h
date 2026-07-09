#pragma once

#include "container/Container.h"
#include "Kernel.h"

class ForceComputation{
    protected:

    struct SPHSolver{
        SPHSolver(Real r):density{r},pressure{r},viscosity{r}{

        }
        DensityKernel density;
        PressureKernel pressure;
        ViscosityKernel viscosity;
    } kernels;

    ScalarField pressure_field;
    Real h2;
    Real constant;
    public:
    ForceComputation(Real r,size_t num_parts);
    void ComputeForces(ParticleContainer& container);
    
    public:
    void ApproximateDensity(ParticleContainer& container);
    void CalculatePressure(ParticleContainer& particles);
};