#include<thrust/scan.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/device_ptr.h>

#include "cudaLaunchers.cuh"
#include "cuNeighborLists.cuh"
#include "order/cuRowMajorOrder.cuh"
#include "order/cuMortonOrder.cuh"

/**
 * KERNEL DECLARATIONS
 * 
 */

__global__
void cuSetLastOffsetElement(int N, int* offset_per_particle, const int* neighbors_per_particle);


void cuCalculateCellOffsets(int N, 
                          int* offset_per_particle,
                          int* neighbors_per_particle,
                          int*& neighbor_list
                            );


/*******************************************************************
 * 
 * LAUNCHER AND TEMPLATE CALLBACK FUNCTIONS
 * 
 *******************************************************************/
void cuCallbackNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id, 
                           int num_cells, int nx, int ny, 
                           int* cell_start, int* cell_end,
                           int* neighbors_per_particle, int* offset_per_particle,
                           int*& neighbor_list,
                           std::string& order
                          ){
    //check for type
    if(order=="morton"){

    cuCountNeighborsPerParticle<cuMortonOrder><<<blocks,threads>>>(N,rx,ry,
                                                    particle_cell_id,
                                                    num_cells,nx,ny,
                                                    cell_start,cell_end,
                                                    neighbors_per_particle);
    cudaDeviceSynchronize();    

    cuCalculateCellOffsets(N,offset_per_particle,neighbors_per_particle, neighbor_list);

    cuAssignIndices<cuMortonOrder><<<blocks,threads>>>(N,
                                                       rx,ry,
                                                       particle_cell_id,
                                                       num_cells,nx,ny,
                                                       cell_start,cell_end,
                                                       offset_per_particle,neighbor_list,neighbors_per_particle);

    }
    if(order=="rowmajor"){
    cuCountNeighborsPerParticle<cuRowMajorOrder><<<blocks,threads>>>(N,rx,ry,
                                                    particle_cell_id,
                                                    num_cells,nx,ny,
                                                    cell_start,cell_end,
                                                    neighbors_per_particle);
    cudaDeviceSynchronize();    

    cuCalculateCellOffsets(N,offset_per_particle,neighbors_per_particle, neighbor_list);

    cuAssignIndices<cuRowMajorOrder><<<blocks,threads>>>(N,
                                                        rx,ry,
                                                        particle_cell_id,
                                                        num_cells,nx,ny,
                                                        cell_start,cell_end,
                                                        offset_per_particle,neighbor_list,neighbors_per_particle);
    
    }
    cudaDeviceSynchronize();

}




/**********************************************************************
 * 
 * KERNEL DEFINITIONS
 * 
 *********************************************************************/


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
