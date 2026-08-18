#pragma once

#include "container/SoAContainer.h"
#include "cuParticles.cuh"
#include "Domain.h"
#include<cuda_runtime.h>

/**
 * Only code that will be called from HOST
 */

void cuDeviceProperties();
void cudaSetConstants(cuConstants& cucst);

/**
 * \brief Copies host to device
 */
void cuCreateContainer(SoAContainer& container, cuParticles& particles);
void cuReadContainer(SoAContainer& container, cuParticles& particles);
void cuWriteContainer(SoAContainer& container, cuParticles& particles);

void cuSetupNeighborLists(SoAContainer& container, cuNeighborList& neighbor_list);

/**
 * TEST
 */

void cuInvertBuffer(int N, 
                    int blocks, int threads,
                    Real * rho, Real* inv_rho);

/**
 * Functions for FindNeighbors
 */
void cuCreateCellBuffers(int num_cells, cuParticles& particles);

void cuKernelFindNeighbors(int blocks, int threads, 
                        int N, int* particle_cell_index,
                        int nx, int ny, Real cell_size,
                        int* cell_start, int* cell_end,
                        Real* rx, Real* ry,
                        Real* vx, Real* vy,
                        Real* fx, Real* fy);

/**
 * Functions for ForceComputation with Sorted + Cell Ranges
 */
void cuKernelComputeForcesSorted(int blocks, int threads, 
                                 int N, 
                                 Real* rx, Real* ry, 
                                 Real* fx, Real* fy, 
                                 Real* vx, Real* vy, 
                                 int* particle_cell_id, 
                                 int num_cells, int nx, int ny, 
                                 int* start, int* end, 
                                 Real* rho, Real* pressure, 
                                 cuConstants kernel_constants);

void cuKernelApproximateDensitySorted(int blocks, int threads, int N, 
                                      Real* rx, Real* ry, 
                                      int* particle_cell_id, 
                                      int num_cells, int nx, int ny, 
                                      int* start, int* end, Real* rho, 
                                      cuConstants kernel_constants);

/**
 * Neighbor List Force Computation Launchers
 */
void cuKernelComputeForcesWithNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           Real* fx, Real* fy,
                           Real* vx, Real* vy,
                           int* particle_cell_id,
                           int num_cells, int nx, int ny,
                           Real* rho, Real* pressure, 
                           int* offset_per_particle, int* neighbor_list,
                           cuConstants kernel_constants
                           );

void cuKernelApproximateDensityWithNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id,
                           int num_cells, int nx, int ny,
                           Real* rho, 
                           int* offset_per_particle, int* neighbor_list,
                           cuConstants kernel_constants
                           );



void cuKernelComputePressureSorted(int blocks, int threads, 
                                   int N, 
                                   Real* rho, Real* pressure);

void cuKernelNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id, 
                           int num_cells, int nx, int ny, 
                           int* cell_start, int* cell_end,
                           int* neighbors_per_particle, int* offset_per_particle,
                           int*& neighbor_list_indices,
                           cuConstants kernel_constants
                          );

/**
 * Utilities section
 */
template<class T> 
void cuSetupBuffer(size_t N, T*& buffer){
    cudaMalloc(&buffer,N*sizeof(T));
}

template<class T>
void cuFreeBuffer(T* buffer){
    cudaFree(buffer);
}

template<class T>
void cuCopyBufferToDevice(size_t N, T* dst, T* src){
    cudaMemcpy(dst,src,N*sizeof(T),cudaMemcpyHostToDevice);
}

template<class T>
void cuCopyBufferToHost(size_t N, T* dst, T* src){
    cudaMemcpy(dst,src,N*sizeof(T),cudaMemcpyDeviceToHost);
}

template<class T>
void cuSetBufferValue(size_t N, T* buffer, T val){
    cudaMemset(buffer,val,N*sizeof(T));
}