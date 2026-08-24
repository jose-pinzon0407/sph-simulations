#pragma once
#include "cuda_callbacks.cuh"

/**
 * Only code that will be called from HOST
 */
void cuDeviceProperties();



void cuCreateContainer(SoAContainer& container, cuParticles& particles);
void cuReadContainer(SoAContainer& container, cuParticles& particles);
void cuWriteContainer(SoAContainer& container, cuParticles& particles);
void cuSetupNeighborLists(SoAContainer& container, cuNeighborList& neighbor_list);
void cuCreateCellBuffers(int num_cells, cuParticles& particles);

void cuInvertBuffer(int N, 
                    int blocks, int threads,
                    Real * rho, Real* inv_rho);







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


/**
 * Utility Functions
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