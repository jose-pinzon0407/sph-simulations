#pragma once

#include "ConfigVariables.hpp"
#include "cuParticles.cuh"

template<class Order>
__global__
void cuCountNeighborsPerParticle(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* neighbors_per_particle){

    int idx = threadIdx.x+blockIdx.x*blockDim.x;
    int offset_i=0;


    if(idx >= N) return;

    Real rix=rx[idx];
    Real riy=ry[idx];

    int cell_id = particle_cell_id[idx];

    // int cix = cell_id%nx;
    // int ciy = (cell_id/nx)%ny;

    int cix=0;
    int ciy=0;
    int ciz=0;
    Order::MapOneToThree(cix,ciy,ciz,cell_id,nx,ny,0);


    for(int dx=-1;dx<=1;++dx){
        for(int dy=-1;dy<=1;++dy){
            int cjx = cix;
            int cjy = ciy;

            cjx += dx;
            if(cjx < 0 || cjx >= nx) continue;

            cjy += dy;

            if(cjy< 0 || cjy>= ny) continue;


            // int cell_j_id = cjx+cjy*nx;
            int cell_j_id = Order::MapThreeToOne(cjx,cjy,0,nx,ny,0);

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            if(cell_j_start<0 || cell_j_end<0){
                continue;
            }

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;
                if(j == idx) continue;

                if(r2>H*H)continue;
                
                offset_i++;

            }

        }
    }

    neighbors_per_particle[idx]=offset_i;
}


template<class Order>
__global__
void cuAssignIndices(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* offset_per_particle, int* neighbor_list, int* neighbors_per_particle){

    int idx = threadIdx.x+blockIdx.x*blockDim.x;


    if(idx >= N) return;

    int offset_i=offset_per_particle[idx];
    int counter=0;
    Real rix=rx[idx];
    Real riy=ry[idx];

    int cell_id = particle_cell_id[idx];

    // int cix = cell_id%nx;
    // int ciy = (cell_id/nx)%ny;

    int cix=0;
    int ciy=0;
    int ciz=0;
    Order::MapOneToThree(cix,ciy,ciz,cell_id,nx,ny,0);

    for(int dx=-1;dx<=1;++dx){
        for(int dy=-1;dy<=1;++dy){
            int cjx = cix;
            int cjy = ciy;

            cjx += dx;
            if(cjx < 0 || cjx >= nx) continue;

            cjy += dy;

            if(cjy< 0 || cjy>= ny) continue;


            // int cell_j_id = cjx+cjy*nx;
            int cell_j_id = Order::MapThreeToOne(cjx,cjy,0,nx,ny,0);

            


            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            if(cell_j_start<0 || cell_j_end<0){
                continue;
            }

            for(int j=cell_j_start;j<cell_j_end;++j){
                if(counter ==neighbors_per_particle[idx]) break;

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;
                if(j == idx) continue;
                if(r2>H*H)continue;
                
                
                neighbor_list[offset_i] = j;
                offset_i++;
                counter++;


            }

        }
    }


}