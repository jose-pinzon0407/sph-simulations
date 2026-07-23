#include "CUDAKernels.cuh"

__global__
void cuApproximateDensity(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, cuConstants kernel_constants);

__global__ 
void cuComputePressure(int N, const Real* rho, Real* pressure);

__global__
void cuComputeForces(int N, Real* rx, Real* ry, Real* fx, Real* fy, Real* vx, Real* vy, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, Real* pressure, cuConstants kernel_constants);

/***************************
 * LAUNCHER FUNCTIONS
 ***************************/

void cuSetBufferToZero(Real* buffer, int val, int count){
    cudaMemset(buffer,val,count*sizeof(Real));
}

void cuSetUpBuffer(int N, Real*& buffer){
    cudaMalloc(&buffer,N*sizeof(Real));
}

void cuSetUpBuffer(int N, int*& buffer){
    cudaMalloc(&buffer,N*sizeof(int));
}

void cuCopyBufferToDevice(Real* dst, Real* src, size_t count){
    cudaMemcpy(dst,src,count*sizeof(Real),cudaMemcpyHostToDevice);
}

void cuCopyBufferToHost(Real* dst, Real* src, size_t count){
    cudaMemcpy(dst,src,count*sizeof(Real),cudaMemcpyDeviceToHost);
}


void cuCopyBufferToDevice(int* dst, int* src, size_t count){
    cudaMemcpy(dst,src,count*sizeof(int),cudaMemcpyHostToDevice);
}

void cuCopyBufferToDevice(int* dst, uint* src, size_t count){
    cudaMemcpy(dst,src,count*sizeof(int),cudaMemcpyHostToDevice);
}

void cuCopyBufferToHost(int* dst, int* src, size_t count){
    cudaMemcpy(dst,src,count*sizeof(int),cudaMemcpyDeviceToHost);
}

void cuFreeMemory(Real* buffer){
    cudaFree(buffer);
}

void cuFreeMemory(int* buffer){
    cudaFree(buffer);
}

void cuLaunchApproximateDensityKernel(int blocks, int threads, int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, cuConstants kernel_constants){
    cuApproximateDensity<<<blocks,threads>>>(N,rx,ry, particle_cell_id,num_cells, nx, ny, start,end,rho, kernel_constants);
    cudaDeviceSynchronize();
}

void cuLaunchComputePressureKernel(int blocks, int threads, int N, Real* rho, Real* pressure){
    cuComputePressure<<<blocks,threads>>>(N,rho,pressure);
    cudaDeviceSynchronize();
}

void cuLaunchComputeForces(int blocks, int threads, int N, Real* rx, Real* ry, Real* fx, Real* fy, Real* vx, Real* vy, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, Real* pressure, cuConstants kernel_constants){

    cuComputeForces<<<blocks,threads>>>(N,rx,ry,fx,fy,vx,vy,particle_cell_id,num_cells,nx,ny,start,end,rho,pressure, kernel_constants);
    cudaDeviceSynchronize();
}

/**
 * KERNEL DECLARATIONS
 */

__device__ 
Real DensityW(Real r2){
    Real h2=H*H;
    return (h2-r2)*(h2-r2)*(h2-r2);
}

__global__
void cuComputePressure(int N, const Real* rho, Real* pressure){
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if(idx>=N) return;

    pressure[idx] = GAS_CONSTANT*(rho[idx]-RHO_0);

}

__global__ 
void cuApproximateDensity(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, cuConstants kernel_constants){
    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real rho_i=0.0;

    int cell_id = particle_cell_id[idx];
    int cix = cell_id%nx;
    int ciy = (cell_id/nx)%ny;

    for(int dx=-1;dx<=1;++dx){
        for(int dy=-1;dy<=1;++dy){
            int cjx = cix;
            int cjy = ciy;

            cjx += dx;
            cjy += dy;

            int cell_j_id = cjx+cjy*nx;

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;

                if(r2>H*H)continue;

                rho_i += kernel_constants.density_constant*DensityW(r2);

            }

        }
    }

    rho[idx]=rho_i;

}

__global__
void cuComputeForces(int N, Real* rx, Real* ry, Real* fx, Real* fy, Real* vx, Real* vy, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, Real* pressure, cuConstants kernel_constants){

    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    Real vix = vx[idx];
    Real viy = vy[idx];

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real rho_i=rho[idx];
    Real p_rho_i2 = pressure[idx]/(rho_i*rho_i);
    
    fx[idx]=0.0;
    fy[idx]=-10.0;

    //Pressure force
    Real fp_x=0.0;
    Real fp_y=0.0;

    //Viscosity force
    Real fv_x=0.0;
    Real fv_y=0.0;


    int cell_id = particle_cell_id[idx];
    int cix = cell_id%nx;
    int ciy = (cell_id/nx)%ny;

    for(int cx=-1;cx<=1;++cx){
        for(int cy=-1;cy<=1;++cy){
            int cjx = cix;
            int cjy = ciy;

            cjx += cx;
            cjy += cy;

            int cell_j_id = cjx+cjy*nx;

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx + dy*dy;
                
                if(r2>kernel_constants.h2)continue;


                Real r1 = sqrt(r2);
                Real inv_r1 = 1.0/r1;
                if(r1>TOL){
                    Real rho_j=rho[j];
                    //pressure forces
                    Real pressure_multiplier = -0.5*(p_rho_i2+pressure[j]/(rho_j*rho_j));

                    Real h1_r12= (kernel_constants.h1-r1)*(kernel_constants.h1-r1);

                    fp_x += -3.0*kernel_constants.pressure_constant*h1_r12*inv_r1*dx*pressure_multiplier;
                    fp_y += -3.0*kernel_constants.pressure_constant*h1_r12*inv_r1*dy*pressure_multiplier;
                }

                Real dvx = vx[j]-vix;
                Real dvy = vy[j]-viy;

                Real viscosity_multiplier = kernel_constants.viscosity_constant/rho[j]*(kernel_constants.h1-r1);

                fv_x += viscosity_multiplier*dvx;
                fv_y += viscosity_multiplier*dvy;


            }

        }
    }

    fx[idx]+= fp_x + fv_x;
    fy[idx]+= fp_y + fv_y;

}
