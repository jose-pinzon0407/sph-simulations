#pragma once

#include "container/SoAContainer.h"


class SoAForceComputation{

    public:
    SoAForceComputation(Real radius):h1{radius},h2{radius*radius},h6{h2*h2*h2}{
        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        std::cout<<"[SoAForceComputation] density constant="<<density_constant<<std::endl;
        std::cout<<"[SoAForceComputation] pressure constant="<<pressure_constant<<std::endl;
    }

    void ComputeForces(SoAContainer& particles);

    public:
    void ApproximateDensity(SoAContainer& particles);
    void CalculatePressure(SoAContainer& particles);
    
    protected:
    Real h1;
    Real h2;
    Real h6;
    ScalarField pressure_field;

    private:
    Real density_constant;
    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }

    Real pressure_constant;
};