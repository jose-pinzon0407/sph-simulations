#pragma once

#include "container/SoAContainer.h"

template<class NeighborSearchType>
class GenericForceComputation{

    public:
    GenericForceComputation(Real radius, GenericNeighborSearch<NeighborSearchType>& ns);

    void Compute(SoAContainer& particles);
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
    GenericNeighborSearch<NeighborSearchType>& neighbor_search;

    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }

};

#include "GenericForceComputation.tpp"