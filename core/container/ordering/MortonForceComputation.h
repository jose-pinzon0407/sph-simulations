#pragma once

#include "container/SoAContainer.h"
#include "container/ordering/MortonNeighborSearch.h"

class MortonSortedForceComputation{
    public:
    MortonSortedForceComputation(Real radius, MortonNeighborSearch& ns):h1{radius},h2{radius*radius},h6{h2*h2*h2},neighbor_search{ns}{
        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        std::cout<<"[SortedForceComputation] density constant="<<density_constant<<std::endl;
        std::cout<<"[SortedForceComputation] pressure constant="<<pressure_constant<<std::endl;
    }

    void ComputeForces(SoAContainer& particles);
    void ApproximateDensity(SoAContainer& particles);
    void CalculatePressure(SoAContainer& particles);


    protected:
    Real h1;
    Real h2;
    Real h6;
    Real density_constant;
    Real pressure_constant;
    ScalarField pressure_field;
    MortonNeighborSearch& neighbor_search;

    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }
    
};