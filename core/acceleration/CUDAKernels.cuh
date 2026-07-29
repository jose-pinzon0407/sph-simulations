#pragma once


#include "ConfigVariables.hpp"

#include<iostream>

/**
 * Provides function declarations for CUDA and THRUST callbacks, to be used 
 * in the main file.
 * All of this functions should be callable from the main file, so ideally no 
 * cuda code can be found here.
 * 
 */

struct cuForceComputation{

};

struct cuConstants{

    cuConstants(Real radius, Real visc):h1{radius},h2{radius*radius},h6{h2*h2*h2},mu{visc}{
        density_constant = 315.0/(64.0*M_PI*std::pow(h1,9));
        pressure_constant = 15.0/(M_PI*std::pow(h1,6.0));
        viscosity_constant = 45.0*mu/(M_PI*h6);
        std::cout<<"[cuKernelConstants] density constant="<<density_constant<<std::endl;
        std::cout<<"[cuKernelConstants] pressure constant="<<pressure_constant<<std::endl;
        std::cout<<"[cuKernelConstants] viscosity constant="<<viscosity_constant<<std::endl;
    }

    Real h1;
    Real h2;
    Real h6;
    Real mu;

    Real density_constant;
    Real pressure_constant;
    Real viscosity_constant;

};  

void cuOutputDeviceProperties();

void cuFindNeighbors(int blocks, int threads
                    ,int N
                    ,int* particle_cell_index
                    ,int* cell_start, int* cell_end
                    ,int nx, int ny, int cell_size
                    ,Real* rx, Real* ry 
                    ,Real* fx, Real* fy
                    ,Real* vx, Real* vy);

void cuSetBufferToZero(Real* buffer, int val, int count);

void cuSetUpBuffer(int N, Real*& buffer);
void cuSetUpBuffer(int N, int*& buffer);

void cuCopyBufferToDevice(Real* dst, Real* src, size_t count);
void cuCopyBufferToHost(Real* dst, Real* src, size_t count);

void cuCopyBufferToDevice(int* dst, int* src, size_t count);
void cuCopyBufferToDevice(int* dst, uint* src, size_t count);
void cuCopyBufferToHost(int* dst, int* src, size_t count);

void cuFreeMemory(Real* buffer);
void cuFreeMemory(int* buffer);

void cuLaunchApproximateDensityKernel(int blocks, int threads, int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, cuConstants kernel_constants);

void cuLaunchComputePressureKernel(int blocks, int threads, int N, Real* rho, Real* pressure);

void cuLaunchComputeForces(int blocks, int threads, int N, Real* rx, Real* ry, Real* fx, Real* fy, Real* vx, Real* vy, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, Real* pressure, cuConstants kernel_constants);

void cuMeasureNeighborsPerParticle(int blocks, int threads, int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, cuConstants kernel_constants);

struct cuBuffers{
    Real* dev_rx=nullptr;
    Real* dev_ry=nullptr;
    Real* dev_rz=nullptr;

    Real* dev_vx=nullptr;
    Real* dev_vy=nullptr;
    Real* dev_vz=nullptr;

    Real* dev_fx=nullptr;
    Real* dev_fy=nullptr;
    Real* dev_fz=nullptr;

    int* cell_start=nullptr;
    int* cell_end=nullptr;
    int* particle_cell_id=nullptr;

    Real* density=nullptr;
    Real* pressure=nullptr;


};