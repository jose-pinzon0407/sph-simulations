#pragma once

#include "ConfigVariables.hpp"


template<class Order>
__global__ 
void cuSetParticlesCellIndex(int N, Real* rx, Real* ry, Real cell_size, int* particle_cell_index, int nx, int ny){

    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    
    if(idx>=N) return;

    int ci0 = floor(rx[idx]/cell_size);
    int ci1 = floor(ry[idx]/cell_size);

    ci0 = min(ci0,nx-1);
    ci1 = min(ci1,ny-1);

    int cell_idx = Order::MapThreeToOne(ci0,ci1,0,nx,ny,0);

    particle_cell_index[idx]=cell_idx; 
}
