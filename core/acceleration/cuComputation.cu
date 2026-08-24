#include "cudaLaunchers.cuh"
#include "cuComputation.cuh"
#include "order/cuRowMajorOrder.cuh"
#include "order/cuMortonOrder.cuh"

/*******************************************************************
 * 
 * KERNEL DECLARATIONS
 * 
 *******************************************************************/

__global__ 
void cuInvertDensity(int N, Real* rho, Real* inv_rho);


__global__ 
void cuComputePressure(int N, const Real* rho, Real* pressure);



/*******************************************************************
 * 
 * LAUNCHER AND TEMPLATE CALLBACK FUNCTIONS
 * 
 *******************************************************************/

void cuInvertBuffer(int N, 
                    int blocks, int threads,
                    Real * rho, Real* inv_rho){
    cuInvertDensity<<<blocks,threads>>>(N,rho,inv_rho);
    cudaDeviceSynchronize();
}

void cuCallbackComputeForcesSorted(int blocks, int threads, 
                                 int N, 
                                 Real* rx, Real* ry, 
                                 Real* fx, Real* fy, 
                                 Real* vx, Real* vy, 
                                 int* particle_cell_id, 
                                 int num_cells, int nx, int ny, 
                                 int* start, int* end, 
                                 Real* rho, Real* pressure, 
                                 cuConstants kernel_constants,
                                 std::string& order){
    //define the type
    if(order=="rowmajor"){
        cuComputeForces<cuRowMajorOrder><<<blocks,threads>>>(N,
        rx,ry,
        fx,fy,
        vx,vy,
        particle_cell_id,
        num_cells,nx,ny,
        start,end,
        rho,pressure,
        kernel_constants);
    }
    if(order=="morton"){
        cuComputeForces<cuMortonOrder><<<blocks,threads>>>(N,
        rx,ry,
        fx,fy,
        vx,vy,
        particle_cell_id,
        num_cells,nx,ny,
        start,end,
        rho,pressure,
        kernel_constants);
    }      
    cudaDeviceSynchronize();
}


void cuCallbackApproximateDensitySorted(int blocks, int threads, int N, 
                                      Real* rx, Real* ry, 
                                      int* particle_cell_id, 
                                      int num_cells, int nx, int ny, 
                                      int* start, int* end, Real* rho, 
                                      cuConstants kernel_constants,
                                      std::string& order){
    //define the type
    if(order=="rowmajor"){
        cuApproximateDensity<cuRowMajorOrder><<<blocks,threads>>>(N,rx,ry, particle_cell_id,num_cells, nx, ny, start,end,rho, kernel_constants);
    }
    if(order=="morton"){
        cuApproximateDensity<cuMortonOrder><<<blocks,threads>>>(N,rx,ry, particle_cell_id,num_cells, nx, ny, start,end,rho, kernel_constants);
    }
    cudaDeviceSynchronize();
}


void cuKernelComputePressureSorted(int blocks, int threads, 
                                   int N, 
                                   Real* rho, Real* pressure){
    cuComputePressure<<<blocks,threads>>>(N,rho,pressure);
    cudaDeviceSynchronize();
}

/*******************************************************************
 * 
 * KERNEL DEFINITIONS
 * 
 *******************************************************************/

__global__
void cuComputePressure(int N, const Real* rho, Real* pressure){
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if(idx>=N) return;

    pressure[idx] = GAS_CONSTANT*(rho[idx]-RHO_0);

}


__global__ 
void cuInvertDensity(int N, Real* rho, Real* inv_rho){

    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    inv_rho[idx] = 1.0/rho[idx];

}