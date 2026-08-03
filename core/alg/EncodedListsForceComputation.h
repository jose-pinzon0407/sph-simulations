#pragma once

#include "SoAForceComputation.h"
#include "NeighborLists.h"

class EncodedListsForceComputation{
    public:
    EncodedListsForceComputation(Real radius, NeighborSearch& ns, NeighborList& nl):h1{radius},h2{radius*radius},h6{h2*h2*h2},neighbor_search{ns},neighbor_list{nl}{
        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        std::cout<<"[EncodedListsForceComputation] density constant="<<density_constant<<std::endl;
        std::cout<<"[EncodedListsForceComputation] pressure constant="<<pressure_constant<<std::endl;
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
    NeighborSearch& neighbor_search;
    NeighborList& neighbor_list;

    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }

};