#include<thrust/sort.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/device_ptr.h>
#include<thrust/iterator/zip_iterator.h>


#include "cudaLaunchers.cuh"
#include "cuFindNeighbors.cuh"
#include "order/cuRowMajorOrder.cuh"
#include "order/cuMortonOrder.cuh"

//We need the thrust sort method, don't want to do ourselves
void cuSortContainerByCellID(int N, int* particle_cell_index, Real* rx, Real* ry, Real* vx, Real* vy, Real* fx, Real* fy);

__global__
void cuSetCellRanges(int N, int* cell_start, int* cell_end, int* particle_cell_index);


/*******************************************************************
 * 
 * LAUNCHER FUNCTIONS
 * 
 *******************************************************************/
void cuCreateCellBuffers(int num_cells, cuParticles& particles){
    cuSetupBuffer(num_cells,particles.cell_end);
    cuSetupBuffer(num_cells,particles.cell_start);

    cuSetBufferValue(num_cells,particles.cell_end,-1);
    cuSetBufferValue(num_cells,particles.cell_start,-1);
}

void cuCallbackFindNeighbors(int blocks, int threads, 
                        int N, int* particle_cell_index,
                        int nx, int ny, Real cell_size,
                        int* cell_start, int* cell_end,
                        Real* rx, Real* ry,
                        Real* vx, Real* vy,
                        Real* fx, Real* fy,
                        std::string& order)
{
    int total_cells=nx*ny;

    //check for type
    if(order=="morton"){
        cuSetParticlesCellIndex<cuMortonOrder><<<blocks,threads>>>(N,rx,ry,cell_size,particle_cell_index,nx,ny);
        cudaDeviceSynchronize();

        cuSortContainerByCellID(N,particle_cell_index,rx,ry,vx,vy,fx,fy);
        cudaDeviceSynchronize();
        cudaMemset(cell_start,-1,total_cells*sizeof(int));
        cudaMemset(cell_end,-1,total_cells*sizeof(int));

        cuSetCellRanges<<<blocks,threads>>>(N,cell_start,cell_end,particle_cell_index);
        cudaDeviceSynchronize();
        return;
    }
    
    if(order=="rowmajor"){
        cuSetParticlesCellIndex<cuRowMajorOrder><<<blocks,threads>>>(N,rx,ry,cell_size,particle_cell_index,nx,ny);
        cudaDeviceSynchronize();

        cuSortContainerByCellID(N,particle_cell_index,rx,ry,vx,vy,fx,fy);
        cudaDeviceSynchronize();
        cudaMemset(cell_start,-1,total_cells*sizeof(int));
        cudaMemset(cell_end,-1,total_cells*sizeof(int));

        cuSetCellRanges<<<blocks,threads>>>(N,cell_start,cell_end,particle_cell_index);
        cudaDeviceSynchronize();
        return;
    }
}


/*******************************************************************
 * 
 * KERNEL DEFINITIONS
 * 
 *******************************************************************/


void cuSortContainerByCellID(int N, int* particle_cell_index, Real* rx, Real* ry, Real* vx, Real* vy, Real* fx, Real* fy){
    thrust::device_ptr<int> ptr_particle_cell_index = thrust::device_pointer_cast<int>(particle_cell_index);

    thrust::device_ptr<Real> ptr_rx = thrust::device_pointer_cast<Real>(rx);
    thrust::device_ptr<Real> ptr_ry = thrust::device_pointer_cast<Real>(ry);

    thrust::device_ptr<Real> ptr_vx = thrust::device_pointer_cast<Real>(vx);
    thrust::device_ptr<Real> ptr_vy = thrust::device_pointer_cast<Real>(vy);
    
    thrust::device_ptr<Real> ptr_fx = thrust::device_pointer_cast<Real>(fx);
    thrust::device_ptr<Real> ptr_fy = thrust::device_pointer_cast<Real>(fy);

    thrust::sort_by_key(ptr_particle_cell_index,ptr_particle_cell_index+N,
                        thrust::make_zip_iterator(
                        thrust::make_tuple(ptr_rx,ptr_ry,ptr_vx,ptr_vy,ptr_fx,ptr_fy)
                        )
    );

}


__global__
void cuSetCellRanges(int N, int* cell_start, int* cell_end, int* particle_cell_index){
    
    int idx = threadIdx.x +  blockIdx.x*blockDim.x;

    if(idx>=N) return;

    int idx_cell = particle_cell_index[idx];

    //if in first cell or different w.r.t. previous
    if(idx==0 || idx_cell != particle_cell_index[idx-1]){
        cell_start[idx_cell]=idx;
    }
    
    if(idx==(N-1) || idx_cell != particle_cell_index[idx+1]){
        cell_end[idx_cell] = idx+1;
    }


}