#include<thrust/sort.h>
#include<thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>
#include "CUDAKernels.cuh"

__global__
void cuSortParticlesIntoCells(int N, Real* rx, Real* ry, Real cell_size, int* particle_cell_index, int nx, int ny);

void cuThrustCallBack(int N, int* particle_cell_index, Real* rx, Real* ry 
                                                     ,Real* fx, Real* fy
                                                     ,Real* vx, Real* vy);

__global__
void cuSetCellRanges(int N, int* cell_start, int* cell_end, int* particle_cell_index);

/************************\\
 * LAUNCHER FUNCTIONS
 *************************/

void cuFindNeighbors(int blocks, int threads
                    ,int N
                    ,int* particle_cell_index
                    ,int* cell_start, int* cell_end
                    ,int nx, int ny, int cell_size
                    ,Real* rx, Real* ry 
                    ,Real* fx, Real* fy
                    ,Real* vx, Real* vy)
{
    cuSortParticlesIntoCells<<<blocks,threads>>>(N,rx,ry,cell_size,
                                        particle_cell_index,
                                        nx,ny);
    cudaDeviceSynchronize();
    cuThrustCallBack(N, particle_cell_index
                                        ,rx,ry
                                        ,fx,fy
                                        ,vx,vy);
    
    //Set cell ranges to -1
    cudaMemset(cell_start,-1,nx*ny*sizeof(int));
    cudaMemset(cell_end,-1,nx*ny*sizeof(int));
    cuSetCellRanges<<<blocks,threads>>>(N,cell_start,cell_end,particle_cell_index);
    cudaDeviceSynchronize();

}



void cuThrustCallBack(int N, int* particle_cell_index, Real* rx, Real* ry 
                                                     ,Real* fx, Real* fy
                                                     ,Real* vx, Real* vy){
    
    //lets not use raw pointers directly
    thrust::device_ptr<int> particl_cell_index_ptr = thrust::device_pointer_cast<int>(particle_cell_index);

    thrust::device_ptr<Real> rx_ptr = thrust::device_pointer_cast<Real>(rx);
    thrust::device_ptr<Real> ry_ptr = thrust::device_pointer_cast<Real>(ry);

    thrust::device_ptr<Real> vx_ptr = thrust::device_pointer_cast<Real>(vx);
    thrust::device_ptr<Real> vy_ptr = thrust::device_pointer_cast<Real>(vy);

    thrust::device_ptr<Real> fx_ptr = thrust::device_pointer_cast<Real>(fx);
    thrust::device_ptr<Real> fy_ptr = thrust::device_pointer_cast<Real>(fy);


    //key, zip iterator, tuple collection
    thrust::sort_by_key(particl_cell_index_ptr,particl_cell_index_ptr+N,
                        thrust::make_zip_iterator(
                            thrust::make_tuple(rx_ptr,ry_ptr
                                                ,fx_ptr,fy_ptr
                                                ,vx_ptr,vy_ptr)
                        ));
}


__global__
void cuSortParticlesIntoCells(int N, Real* rx, Real* ry, Real cell_size, int* particle_cell_index, int nx, int ny){

    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    
    if(idx>=N) return;

    int ci0 = floor(rx[idx]/cell_size);
    int ci1 = floor(ry[idx]/cell_size);

    ci0 = min(ci0,nx-1);
    ci1 = min(ci1,ny-1);

    int cell_idx = ci0 + ci1*nx;
    particle_cell_index[idx] = cell_idx;
    

}

__global__
void cuSetCellRanges(int N, int* cell_start, int* cell_end, int* particle_cell_index){

    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if(idx>=N) return;

    int idx_cell = particle_cell_index[idx];

    if(idx ==0 || idx_cell != particle_cell_index[idx-1]){
        cell_start[idx_cell] = idx;
    }

    if(idx == (N-1) || idx_cell != particle_cell_index[idx+1]){
        cell_end[idx_cell] = idx+1;
    }

}