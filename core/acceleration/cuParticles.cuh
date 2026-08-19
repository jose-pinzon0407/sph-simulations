#pragma once

#include "ConfigVariables.hpp"


struct cuParticles{
    Real* dev_rx=nullptr;
    Real* dev_ry=nullptr;
    Real* dev_rz=nullptr;

    Real* dev_vx=nullptr;
    Real* dev_vy=nullptr;
    Real* dev_vz=nullptr;
 
    Real* dev_fx=nullptr;
    Real* dev_fy=nullptr;
    Real* dev_fz=nullptr;

    int* particle_cell_id=nullptr;

    Real* density=nullptr;
    Real* pressure=nullptr;

    //for neighbor lists only
    int* cell_start=nullptr;
    int* cell_end=nullptr;


};

struct cuGridInformation{
    Real cell_size=0;
    int nx, ny, nz;

    int GetTotalCells(){
        return nx*ny;
    }
};

struct cuNeighborList{
    int* dev_neighbors_per_particle;//size N
    int* dev_offset_per_particle;//size N+1
    int* dev_neighbor_indices;// size N
};

struct cuConstants{

    cuConstants(Real radius, Real visc):h1{radius},h2{radius*radius},h6{h2*h2*h2},mu{visc}{
        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        viscosity_constant = 45.0*mu/(M_PI*h6);
        std::cout<<"[cuKernelConstants] density constant="<<density_constant<<std::endl;
        std::cout<<"[cuKernelConstants] pressure constant="<<pressure_constant<<std::endl;
        std::cout<<"[cuKernelConstants] viscosity constant="<<viscosity_constant<<std::endl;

        std::cout<<"[cuKernelConstants] h1="<<h1<<std::endl;
        std::cout<<"[cuKernelConstants] h2="<<h2<<std::endl;
        std::cout<<"[cuKernelConstants] h6="<<h6<<std::endl;

    }

    Real h1;
    Real h2;
    Real h6;
    Real mu;

    Real density_constant;
    Real pressure_constant;
    Real viscosity_constant;

};  