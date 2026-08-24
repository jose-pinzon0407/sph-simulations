#include "cuIntegratorLauncher.cuh"

__global__
void cuUpdatePositions(int N,
                       Real ubx, Real uby, 
                       Real lbx, Real lby, 
                       Real dt,
                       Real* rx, Real* ry, 
                       Real* vx, Real* vy);

__global__
void cuUpdateVelocities(int N,
                        Real dt,
                        Real* fx, Real* fy, 
                        Real* vx, Real* vy);

void cuCallbackIntegrate(int blocks, int threads,
                         int N,
                         Real ubx, Real uby,
                         Real lbx, Real lby,
                         Real dt,
                         Real* fx, Real* fy,
                         Real* rx, Real* ry,
                         Real* vx, Real* vy){
    cuKernelUpdateVelocities(blocks,threads,N,dt,fx,fy,vx,vy);
    cuKernelUpdatePositions(blocks,threads,N,ubx,uby,lbx,lby,dt,rx,ry,vx,vy);
    cudaDeviceSynchronize();
}

void cuKernelUpdatePositions(int blocks, int threads, 
                             int N,
                             Real ubx, Real uby, 
                             Real lbx, Real lby, 
                             Real dt,
                             Real* rx, Real* ry, 
                             Real* vx, Real* vy
                             ){
//
    cuUpdatePositions<<<blocks,threads>>>(N,ubx,uby,lbx,lby,
                                          dt,
                                          rx,ry,
                                          vx,vy  
                                          );
}


void cuKernelUpdateVelocities(int blocks, int threads, 
                             int N,
                             Real dt,
                             Real* fx, Real* fy, 
                             Real* vx, Real* vy
                             ){
    //
    cuUpdateVelocities<<<blocks,threads>>>(N,dt,
                                           fx,fy,
                                           vx,vy
                                           );
}

__global__
void cuUpdatePositions(int N,
                             Real ubx, Real uby, 
                             Real lbx, Real lby, 
                             Real dt,
                             Real* rx, Real* ry, 
                             Real* vx, Real* vy){
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if(idx >= N) return;
    
    Real& rix = rx[idx];
    Real& riy = ry[idx];

    Real& vix = vx[idx];
    Real& viy = vy[idx];

    //update idx position
    rx[idx] += dt*vx[idx];
    ry[idx] += dt*vy[idx];

    //apply boundary conditions

    if(rix>ubx){
        rix=ubx;
        vix *= -1.0;
    }

    if(riy>uby){
        riy=uby;
        viy *= -1.0;
    }


    if(rix<lbx){
        rix=lbx;
        vix *= -1.0;
    }

    if(riy<lby){
        riy=lby;
        viy *= -1.0;
    }


}

__global__
void cuUpdateVelocities(int N,
                        Real dt,
                        Real* fx, Real* fy, 
                        Real* vx, Real* vy){
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if(idx >= N) return;

    //update idx velocity
    vx[idx] += dt*fx[idx];
    vy[idx] += dt*fy[idx];
    
                            
}