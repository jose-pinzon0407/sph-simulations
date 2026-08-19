#include <thrust/scan.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/device_ptr.h>

#include "cudaLaunchers.cuh"



/**
 * KERNEL DECLARATIONS
 * 
 */

__global__
void cuSetLastOffsetElement(int N, int* offset_per_particle, const int* neighbors_per_particle);

__global__
void cuCountNeighborsPerParticle(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* rho);

void cuCalculateCellOffsets(int N, 
                          int* offset_per_particle,
                          int* neighbors_per_particle,
                          int*& neighbor_list
                            );

__global__
void cuAssignIndices(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* offset_per_particle, int* neighbor_list);


/*******************************************************************
 * 
 * LAUNCHER FUNCTIONS
 * 
 *******************************************************************/

void cuKernelNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id, 
                           int num_cells, int nx, int ny, 
                           int* cell_start, int* cell_end,
                           int* neighbors_per_particle, int* offset_per_particle, int*& neighbor_list                          ){

    cuCountNeighborsPerParticle<<<blocks,threads>>>(N,rx,ry,
                                                particle_cell_id,
                                                num_cells,nx,ny,
                                                cell_start,cell_end,
                                                neighbors_per_particle);
    cudaDeviceSynchronize();    

    cuCalculateCellOffsets(N,offset_per_particle,neighbors_per_particle, neighbor_list);

    cuAssignIndices<<<blocks,threads>>>(N,rx,ry,particle_cell_id,num_cells,nx,ny,cell_start,cell_end,offset_per_particle,neighbor_list);
    cudaDeviceSynchronize();
} 


/**********************************************************************
 * 
 * KERNEL DEFINITIONS
 * 
 *********************************************************************/


 __global__
void cuCountNeighborsPerParticle(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* neighbors_per_particle
){

    int idx = threadIdx.x+blockIdx.x*blockDim.x;
    int offset_i=0;


    if(idx >= N) return;

    Real rix=rx[idx];
    Real riy=ry[idx];

    int cell_id = particle_cell_id[idx];

    int cix = cell_id%nx;
    int ciy = (cell_id/nx)%ny;


    for(int dx=-1;dx<=1;++dx){
        for(int dy=-1;dy<=1;++dy){
            int cjx = cix;
            int cjy = ciy;

            cjx += dx;
            if(cjx < 0 || cjx >= nx) continue;

            cjy += dy;

            if(cjy< 0 || cjy>= ny) continue;


            int cell_j_id = cjx+cjy*nx;

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;

                if(r2>H*H)continue;
                
                offset_i++;

            }

        }
    }

    neighbors_per_particle[idx]=offset_i;

}

__global__
void cuSetLastOffsetElement(int N, int* offset_per_particle, const int* neighbors_per_particle){
    if(threadIdx.x==0 && blockIdx.x==0){
        offset_per_particle[N]=offset_per_particle[N-1]+neighbors_per_particle[N-1];
    }
}

void cuCalculateCellOffsets(int N, 
                          int* offset_per_particle,
                          int* neighbors_per_particle,
                          int*& neighbor_list
                            ){

    
    
    cuFreeBuffer(neighbor_list);
    thrust::device_ptr<int> ptr_opp = thrust::device_pointer_cast<int>(offset_per_particle);
    thrust::device_ptr<int> ptr_npp = thrust::device_pointer_cast<int>(neighbors_per_particle);

    //perform thrust exclusive scan (parallel scan)
    thrust::exclusive_scan(ptr_npp,ptr_npp+N,ptr_opp);

    int num_neighbors = thrust::reduce( ptr_npp,ptr_npp+N,0,thrust::plus<int>());

    cuSetLastOffsetElement<<<1,1>>>(N,offset_per_particle,neighbors_per_particle);
    cudaDeviceSynchronize();
    // std::cout<<"Total neighbors= "<<num_neighbors<<std::endl;
    cuSetupBuffer(num_neighbors, neighbor_list);
    
    
}

__global__
void cuAssignIndices(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, int* offset_per_particle, int* neighbor_list){

    int idx = threadIdx.x+blockIdx.x*blockDim.x;


    if(idx >= N) return;

    int offset_i=offset_per_particle[idx];

    Real rix=rx[idx];
    Real riy=ry[idx];

    int cell_id = particle_cell_id[idx];

    int cix = cell_id%nx;
    int ciy = (cell_id/nx)%ny;


    for(int dx=-1;dx<=1;++dx){
        for(int dy=-1;dy<=1;++dy){
            int cjx = cix;
            int cjy = ciy;

            cjx += dx;
            if(cjx < 0 || cjx >= nx) continue;

            cjy += dy;

            if(cjy< 0 || cjy>= ny) continue;


            int cell_j_id = cjx+cjy*nx;

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;

                if(r2>H*H)continue;
                
                neighbor_list[offset_i] = j;
                offset_i++;

            }

        }
    }


}