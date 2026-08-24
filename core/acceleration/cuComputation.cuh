#pragma once

#include "ConfigVariables.hpp"
#include "cuParticles.cuh"

template<class Order>
__global__
void cuApproximateDensity(int N, Real* rx, Real* ry, 
                          int* particle_cell_id, 
                          int num_cells, int nx, int ny, 
                          int* start, int* end, Real* rho, 
                          cuConstants kernel_constants){
    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    Real rix = rx[idx];
    Real riy = ry[idx];

    Real rho_i=0.0;

    int cell_id = particle_cell_id[idx];

    int cix=0;
    int ciy=0;
    int ciz=0;
    Order::MapOneToThree(cix,ciy,ciz,cell_id,nx,ny,0);

    // int cix = cell_id%nx;
    // int ciy = (cell_id/nx)%ny;

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

            for(int j=cell_j_start;j<cell_j_end;++j){

                Real dx = rix-rx[j];
                Real dy = riy-ry[j];

                Real r2 = dx*dx +dy*dy;

                if(r2>H*H)continue;

                rho_i += kernel_constants.density_constant*(kernel_constants.h2-r2)*(kernel_constants.h2-r2)*(kernel_constants.h2-r2);

            }

        }
    }

    rho[idx]=rho_i;
}

template<class Order>
__global__
void cuComputeForces(int N, 
                    Real* rx, Real* ry, 
                    Real* fx, Real* fy, 
                    Real* vx, Real* vy, 
                    int* particle_cell_id, 
                    int num_cells, int nx, int ny, 
                    int* start, int* end, 
                    Real* rho, Real* pressure,
                    cuConstants constants
                    ){


    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;

    Real vix = vx[idx];
    Real viy = vy[idx];

    Real rix = rx[idx];
    Real riy = ry[idx];

    // Real p_rho_i2 = pressure[idx]/(rho[idx]*rho[idx]);
    Real p_rho_i2 = pressure[idx]*(rho[idx]*rho[idx]);
    
    fx[idx]=0.0;
    fy[idx]=-10.0;

    //Pressure force
    Real fp_x=0.0;
    Real fp_y=0.0;

    //Viscosity force
    Real fv_x=0.0;
    Real fv_y=0.0;


    int cell_id = particle_cell_id[idx];
    int cix=0;
    int ciy=0;
    int ciz=0;
    Order::MapOneToThree(cix,ciy,ciz,cell_id,nx,ny,0);

    for(int cx=-1;cx<=1;++cx){
        for(int cy=-1;cy<=1;++cy){

            int cjx = cix;
            int cjy = ciy;

            cjx += cx;
            if(cjx < 0 || cjx >= nx) continue;

            cjy += cy;

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

                Real r2 = dx*dx + dy*dy;
                              
                if(r2>constants.h2 || idx==j)continue;


                Real pj= pressure[j];  
                // Real inv_rho_j = 1.0/rho[j];
                Real inv_rho_j = rho[j];



                {
                float r2_f =(float)r2;
                float inv_r1_f1 = rsqrtf(r2_f);
                float inv_r1 = inv_r1_f1;
                // Real inv_r1 = rsqrt(r2);
                Real r1 = r2*inv_r1;

                Real pressure_multiplier = -0.5*(p_rho_i2+pj*(inv_rho_j*inv_rho_j));

                
                pressure_multiplier *= -3.0*(constants.h1-r1)*(constants.h1-r1)*constants.pressure_constant*inv_r1;

 

                fp_x = fma(pressure_multiplier,dx,fp_x);
                fp_y = fma(pressure_multiplier,dy,fp_y);


                

                Real dvx = vx[j]-vix;
                Real dvy = vy[j]-viy;

                Real viscosity_multiplier = constants.viscosity_constant*inv_rho_j*(constants.h1-r1);



                fv_x = fma(viscosity_multiplier,dvx,fv_x);
                fv_y = fma(viscosity_multiplier,dvy,fv_y);


                }



            }

        }
    }

    fx[idx]+= fp_x + fv_x;
    fy[idx]+= fp_y + fv_y;

}