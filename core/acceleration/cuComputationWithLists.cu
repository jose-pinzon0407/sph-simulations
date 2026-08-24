#include "cudaLaunchers.cuh"
/*******************************************************************
 * 
 * KERNEL DECLARATION
 * 
 *******************************************************************/

__global__
void cuApproximateDensityWithLists(int N, 
                    Real* rx, Real* ry, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny,  
                    Real* rho, 
                    int* offset_per_particle, int* neighbor_list,
                    cuConstants kernel_constants);

__global__
void cuComputeForcesWithLists(int N, 
                    Real* rx, Real* ry, 
                    Real* fx, Real* fy, 
                    Real* vx, Real* vy, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny,  
                    Real* rho, Real* pressure, 
                    int* offset_per_particle, int* neighbor_list,
                    cuConstants kernel_constants);

/*******************************************************************
 * 
 * LAUNCHER FUNCTIONS
 * 
 *******************************************************************/


void cuKernelApproximateDensityWithNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id,
                           int num_cells, int nx, int ny,
                           Real* rho, 
                           int* offset_per_particle, int* neighbor_list,
                           cuConstants kernel_constants
){

    cuApproximateDensityWithLists<<<blocks,threads>>>(N,rx,ry,
                                                      particle_cell_id,
                                                      num_cells,nx,ny,
                                                      rho,
                                                      offset_per_particle,neighbor_list,
                                                      kernel_constants);
    cudaDeviceSynchronize();

}

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
){
    cuComputeForcesWithLists<<<blocks,threads>>>(N,
                            rx,ry,
                            fx,fy,
                            vx,vy,
                            particle_cell_id,
                            num_cells,nx,ny,
                            rho,pressure,
                            offset_per_particle,neighbor_list,
                            kernel_constants
                        );
    cudaDeviceSynchronize();
}


/*******************************************************************
 * 
 * KERNEL DEFINITIONS
 * 
 *******************************************************************/

__global__
void cuApproximateDensityWithLists(int N, 
                    Real* rx, Real* ry, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny,  
                    Real* rho, 
                    int* offset_per_particle, int* neighbor_list,
                    cuConstants kernel_constants){
    //Start density approximation
    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    int offset_i = offset_per_particle[idx];
    int end_i = offset_per_particle[idx+1];

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real rho_i=0.0;



    for(int j=offset_i;j<end_i;++j){
        int idx_j=neighbor_list[j];

        Real dx = rix-rx[idx_j];
        Real dy = riy-ry[idx_j];

        Real r2 = dx*dx + dy*dy;


        rho_i += kernel_constants.density_constant*(kernel_constants.h2-r2)*(kernel_constants.h2-r2)*(kernel_constants.h2-r2);

    }

    rho_i += kernel_constants.density_constant*(kernel_constants.h2)*(kernel_constants.h2)*(kernel_constants.h2);


    rho[idx]=rho_i;


}

__global__
void cuComputeForcesWithLists(int N, 
                    Real* rx, Real* ry, 
                    Real* fx, Real* fy, 
                    Real* vx, Real* vy, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny,  
                    Real* rho, Real* pressure, 
                    int* offset_per_particle, int* neighbor_list,
                    cuConstants kernel_constants){
    //Start force computation
    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    int offset_i = offset_per_particle[idx];
    int end_i = offset_per_particle[idx+1];
    
    Real vix = vx[idx];
    Real viy = vy[idx];

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real rho_i=rho[idx];
    Real p_rho_i2 = pressure[idx]*(rho_i*rho_i);
    
    fx[idx]=0.0;
    fy[idx]=-10.0;

    //Pressure force
    Real fp_x=0.0;
    Real fp_y=0.0;

    //Viscosity force
    Real fv_x=0.0;
    Real fv_y=0.0;

    for(int j=offset_i;j<end_i;++j){
        int idx_j=neighbor_list[j];
        // if(idx_j == idx) continue;
        Real dx = rix-rx[idx_j];
        Real dy = riy-ry[idx_j];

        Real r2 = dx*dx + dy*dy;
        Real pj= pressure[idx_j];  
        Real rho_j=rho[idx_j];

        float r2_f =(float)r2;
        float inv_r1_f1 = rsqrtf(r2_f);
        float inv_r1 = inv_r1_f1;
        // Real inv_r1 = rsqrt(r2);
        Real r1 = r2*inv_r1;

        

        {


            Real pressure_multiplier = 1.5*(kernel_constants.h1-r1)*(kernel_constants.h1-r1)*kernel_constants.pressure_constant*inv_r1*(p_rho_i2+pj*(rho_j*rho_j));
            

            fp_x = fma(pressure_multiplier,dx,fp_x);
            fp_y = fma(pressure_multiplier,dy,fp_y);

                

            Real dvx = vx[idx_j]-vix;
            Real dvy = vy[idx_j]-viy;

            Real viscosity_multiplier = kernel_constants.viscosity_constant/rho_j*(kernel_constants.h1-r1);

            fv_x = fma(viscosity_multiplier,dvx,fv_x);
            fv_x = fma(viscosity_multiplier,dvy,fv_y);


        }

    }

    fx[idx]+= fp_x + fv_x;
    fy[idx]+= fp_y + fv_y;

}