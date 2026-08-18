#include "cudaLaunchers.cuh"

/**
 * CONSTANTS
 */
__constant__ Real density_constant;
__constant__ Real pressure_constant;
__constant__ Real viscosity_constant;
__constant__ Real h1;
__constant__ Real h2;


void cudaSetConstants(cuConstants& cucst){
        cudaMemcpyToSymbol(pressure_constant,&cucst.pressure_constant,sizeof(Real));
        cudaMemcpyToSymbol(viscosity_constant,&cucst.viscosity_constant,sizeof(Real));
        cudaMemcpyToSymbol(density_constant,&cucst.density_constant,sizeof(Real));
        cudaMemcpyToSymbol(h1,&cucst.h1,sizeof(Real));
        cudaMemcpyToSymbol(h2,&cucst.h2,sizeof(Real));
}


/*******************************************************************
 * 
 * KERNEL DECLARATIONS
 * 
 *******************************************************************/

__global__ 
void cuInvertDensity(int N, Real* rho, Real* inv_rho);

__global__
void cuApproximateDensity(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, cuConstants kernel_constants);

__global__ 
void cuComputePressure(int N, const Real* rho, Real* pressure);

__global__
void cuComputeForces(int N, 
                    Real* rx, Real* ry, 
                    Real* fx, Real* fy, 
                    Real* vx, Real* vy, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny, 
                    int* start, int* end, 
                    Real* rho, Real* pressure);

__device__
Real DensityW(Real r2);

/*******************************************************************
 * 
 * LAUNCHER FUNCTIONS
 * 
 *******************************************************************/

void cuInvertBuffer(int N, 
                    int blocks, int threads,
                    Real * rho, Real* inv_rho){
    cuInvertDensity<<<blocks,threads>>>(N,rho,inv_rho);
}

void cuKernelComputeForcesSorted(int blocks, int threads, int N, Real* rx, Real* ry, Real* fx, Real* fy, Real* vx, Real* vy, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, Real* rho, Real* pressure, cuConstants kernel_constants){
    cuComputeForces<<<blocks,threads>>>(N,
        rx,ry,
        fx,fy,
        vx,vy,
        particle_cell_id,
        num_cells,nx,ny,
        start,end,
        rho,pressure);
    cudaDeviceSynchronize();
}

void cuKernelApproximateDensitySorted(int blocks, int threads, int N, 
                                      Real* rx, Real* ry, 
                                      int* particle_cell_id, 
                                      int num_cells, int nx, int ny, 
                                      int* start, int* end, 
                                      Real* rho, 
                                      cuConstants kernel_constants){

    //call kernel
    cuApproximateDensity<<<blocks,threads>>>(N,rx,ry, particle_cell_id,num_cells, nx, ny, start,end,rho, kernel_constants);
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

                rho_i += kernel_constants.density_constant*DensityW(r2);

            }

        }
    }

    rho[idx]=rho_i;

}


__global__
void cuComputeForces(int N, 
    Real* rx, Real* ry, 
    Real* fx, Real* fy, 
    Real* vx, Real* vy, 
    int* particle_cell_id, 
    int num_cells, int nx, int ny, 
    int* start, int* end, 
    Real* rho, Real* pressure
    ){

    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    Real vix = vx[idx];
    Real viy = vy[idx];

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real p_rho_i2 = pressure[idx]/(rho[idx]*rho[idx]);
    
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
            if(cjx < 0 || cjx >= nx) continue;

            cjy += cy;

            if(cjy< 0 || cjy>= ny) continue;

            int cell_j_id = cjx+cjy*nx;

            int cell_j_start = start[cell_j_id];
            int cell_j_end = end[cell_j_id];

            if(cell_j_start<0 || cell_j_end<0){
                continue;
            }

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx + dy*dy;
                              
                if(r2>h2 || idx==j)continue;


                Real pj= pressure[j];  
                Real inv_rho_j = 1.0/rho[j];


                // if(r2<=TOL*TOL) continue;

                {
                
                Real inv_r1 = rsqrt(r2);
                Real r1 = r2*inv_r1;

                // Real pressure_multiplier = -0.5*(p_rho_i2+pj/(rho_j*rho_j));
                Real pressure_multiplier = -0.5*(p_rho_i2+pj*(inv_rho_j*inv_rho_j));

                // Real pressure_multiplier = -0.5*(pressure[idx]/(rho[idx]*rho[idx])+pressure[j]/(rho[j]*rho[j]));

                // Real h1_r12= (kernel_constants.h1-r1)*(kernel_constants.h1-r1);
                
                pressure_multiplier *= -3.0*(h1-r1)*(h1-r1)*pressure_constant*inv_r1;

                // fx[idx] += -3.0*kernel_constants.pressure_constant*h1_r12*inv_r1*dx*pressure_multiplier;
                // fy[idx] += -3.0*kernel_constants.pressure_constant*h1_r12*inv_r1*dy*pressure_multiplier;
                
                // fx[idx] += pressure_multiplier*dx;
                // fy[idx] += pressure_multiplier*dy;

                // fp_x += pressure_multiplier*dx;
                // fp_y += pressure_multiplier*dy;

                fp_x = fma(pressure_multiplier,dx,fp_x);
                fp_y = fma(pressure_multiplier,dy,fp_y);


                
                // fx[idx] += -3.0*kernel_constants.pressure_constant*(kernel_constants.h1-r1)*(kernel_constants.h1-r1)*inv_r1*(-0.5*(pressure[idx]/(rho[idx]*rho[idx])+pressure[j]/(rho[j]*rho[j])))*dx;
                // fx[idx] += -3.0*kernel_constants.pressure_constant*(kernel_constants.h1-r1)*(kernel_constants.h1-r1)*inv_r1*(-0.5*(pressure[idx]/(rho[idx]*rho[idx])+pressure[j]/(rho[j]*rho[j])))*dy;
                

                Real dvx = vx[j]-vix;
                Real dvy = vy[j]-viy;

                // Real viscosity_multiplier = kernel_constants.viscosity_constant/rho_j*(kernel_constants.h1-r1);
                Real viscosity_multiplier = viscosity_constant*inv_rho_j*(h1-r1);


                // fx[idx] += viscosity_multiplier*dvx;
                // fy[idx] += viscosity_multiplier*dvy;

                // fv_x += viscosity_multiplier*dvx;
                // fv_y += viscosity_multiplier*dvy;

                fv_x = fma(viscosity_multiplier,dvx,fv_x);
                fv_y = fma(viscosity_multiplier,dvy,fv_y);



                // fx[idx] += kernel_constants.viscosity_constant/rho[j]*(kernel_constants.h1-r1)*(vx[j]-vx[idx]);
                // fx[idx] += kernel_constants.viscosity_constant/rho[j]*(kernel_constants.h1-r1)*(vy[j]-vy[idx]);

                }

                // fv_x += dvx;
                // fv_y += dvy;



            }

        }
    }

    fx[idx]+= fp_x + fv_x;
    fy[idx]+= fp_y + fv_y;

}



__global__ 
void cuInvertDensity(int N, Real* rho, Real* inv_rho){

    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    inv_rho[idx] = 1.0/rho[idx];

}