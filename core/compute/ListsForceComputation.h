#pragma once

#include "container/SoAContainer.h"

class ListsForceComputation{

    public:
    ListsForceComputation(Real radius, std::vector<std::vector<int>>& nl):h1{radius},h2{radius*radius},h6{h2*h2*h2},neighbor_list{nl}{

        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        std::cout<<"[ListsForceComputation] density constant="<<density_constant<<std::endl;
        std::cout<<"[ListsForceComputation] pressure constant="<<pressure_constant<<std::endl;
    }
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
    std::vector<std::vector<int>>& neighbor_list;

    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }

};


class MPListsForceComputation{

    public:
    MPListsForceComputation(Real radius, std::vector<int>& nl, std::vector<int>& opp):h1{radius},h2{radius*radius},h6{h2*h2*h2},neighbor_list{nl},offset_per_particle{opp}{

        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        std::cout<<"[MPListsForceComputation] density constant="<<density_constant<<std::endl;
        std::cout<<"[MPListsForceComputation] pressure constant="<<pressure_constant<<std::endl;
    }
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
    std::vector<int>& neighbor_list;
    std::vector<int>& offset_per_particle;

    Real DensityW(Real r2){
        return density_constant*(h2-r2)*(h2-r2)*(h2-r2);
    }

};