#include "Kernel.h"
#include "Domain.h"

#include "container/MeasuringFunctions.h"
#include "container/Lattice.h"
#include "container/Container.h"
#include "container/SoAContainer.h"
#include "container/SortedContainer.h"

#include "alg/ForceComputation.h"
#include "alg/SoAForceComputation.h"
#include "alg/NeighborSearch.h"

#include "integrator/Integrator.h"


#include "utils/CheckPointWriter.h"
#include "utils/Timer.h"
#include "acceleration/CUDAKernels.cuh"

// #ifdef GPU
// #include "acceleration/CUDAKernels.cuh"
// #endif

std::string StepHeader(int step, Container& particles, Domain& domain);

Timer force_computation{"force_computation"};
Timer density_approximation{"density_approximation"};
Timer container_build{"container_build"};
Timer main_loop{"main_loop"};
Timer cu_copy_buffer{"cu_copy"};
Timer integrate_operations{"integrate_operations"};
Timer checkpointing{"checkpointing"};

int main(int argc, char** argv){
    
    #ifdef GPU
    std::cout<<"[Build] Using GPU acceleration\n";
    cuOutputDeviceProperties();
    #endif
    
    int print_frequency=10;
    int steps=10;
    Real dt=0.005;
    Domain domain;
    domain.h_corner=Vec3d {40,50,0};
    domain.l_corner=Vec3d {0,0,0};
    NeighborSearch neighbor_search{domain,H};

    std::cout<<"[Domain] Length:"<<domain.GetLength()<<std::endl;

    // ParticleContainer container;
    // SoAContainer container;
    SortedSoAContainer container;
    ViscosityKernel kernel{H};
    
    Integrator integrator{dt,&domain};

    LatticeCreator lattice;
    lattice.yVol=0.9;
    lattice.xVol=.9;
    lattice.dX=0.125;
    lattice.dY=0.125;

    lattice.CreateLattice(domain,container);
    neighbor_search.FindNeighbors(container);

    CheckPointWriter checkpoint{"positions_step_",0};
    checkpoint.WriteCheckPoint(0,container);

    //Device memory pointers
    cuBuffers dev_buffers;
    int N=container.Size();

    //positions
    cuSetUpBuffer(N,dev_buffers.dev_rx);
    cuSetUpBuffer(N,dev_buffers.dev_ry);
    //velocities
    cuSetUpBuffer(N,dev_buffers.dev_vx);
    cuSetUpBuffer(N,dev_buffers.dev_vy);
    //forces
    cuSetUpBuffer(N,dev_buffers.dev_fx);
    cuSetUpBuffer(N,dev_buffers.dev_fy);
    cuSetUpBuffer(N,dev_buffers.particle_cell_id);

    int num_cells = neighbor_search.GetTotalCells();
    int nx = neighbor_search.GetCellsInDirectionX();
    int ny = neighbor_search.GetCellsInDirectionY();

    cuSetUpBuffer(num_cells,dev_buffers.cell_start);
    cuSetUpBuffer(num_cells,dev_buffers.cell_end);

    cuSetUpBuffer(N,dev_buffers.density);
    cuSetUpBuffer(N,dev_buffers.pressure);

    //auxiliary
    std::vector<Real> host_density;
    host_density.resize(N);
    std::fill(host_density.begin(),host_density.end(),0.0);

    std::vector<Real> host_pressure;
    host_pressure.resize(N);
    std::fill(host_pressure.begin(),host_pressure.end(),0.0);
    //Device memory pointer

    //Create objects for CUDA
    cuConstants cuda_constants{H,MU};


    //Launch kernel configuration
    int threads = 256;
    int blocks = (N+threads-1)/threads;
    std::cout<<"[GPU] Launching: "<<blocks<<" blocks with "<<threads<<" threads grid\n";
    std::cout<<"[GPU] Total threads: "<<blocks*threads<<std::endl;


    // SoAForceComputation computation{H};
    SortedForceComputation computation{H,neighbor_search};

    Particles& particles = container.GetParticles();
    main_loop.Start();
    for(int i=1;i<steps;++i){
        
        // std::cout<<StepHeader(i,container,domain);
        integrate_operations.Start();
        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);
        integrate_operations.Stop();

        checkpointing.Start();
        if(i%print_frequency==0){
            checkpoint.WriteCheckPoint(i,container);
        }
        checkpointing.Stop();

        //I. the particle positions
        cu_copy_buffer.Start();
        cuCopyBufferToDevice(dev_buffers.dev_rx,particles.rx.data(),N);
        cuCopyBufferToDevice(dev_buffers.dev_ry,particles.ry.data(),N);
        // cuCopyBufferToDevice(dev_buffers.particle_cell_id,particles.cell_id.data(),N);
        cu_copy_buffer.Stop();

        // copy velocities
        cuCopyBufferToDevice(dev_buffers.dev_vx,particles.vx.data(),N);
        cuCopyBufferToDevice(dev_buffers.dev_vy,particles.vy.data(),N);

        container_build.Start();
        cuFindNeighbors(blocks,threads,N,dev_buffers.particle_cell_id,dev_buffers.cell_start,dev_buffers.cell_end,nx,ny,cuda_constants.h1,dev_buffers.dev_rx,dev_buffers.dev_ry,dev_buffers.dev_fx,dev_buffers.dev_fy,dev_buffers.dev_vx,dev_buffers.dev_vy);
        container_build.Stop();

        //0.prepare cuda copies
        cuSetBufferToZero(dev_buffers.density,0.0,N);
        cuSetBufferToZero(dev_buffers.pressure,0.0,N);


        density_approximation.Start();
        cuLaunchApproximateDensityKernel(blocks,threads,N,
                                        dev_buffers.dev_rx,dev_buffers.dev_ry, dev_buffers.particle_cell_id,
                                        num_cells,nx,ny,
                                        dev_buffers.cell_start,dev_buffers.cell_end,dev_buffers.density,
                                        cuda_constants
                                        );

        cuLaunchComputePressureKernel(blocks,threads,N,dev_buffers.density,dev_buffers.pressure);
        density_approximation.Stop();


        force_computation.Start();
        cuLaunchComputeForces(blocks,threads,N,
                              dev_buffers.dev_rx,dev_buffers.dev_ry,
                              dev_buffers.dev_fx,dev_buffers.dev_fy,
                              dev_buffers.dev_vx,dev_buffers.dev_vy,dev_buffers.particle_cell_id,
                              num_cells,nx,ny,
                              dev_buffers.cell_start,dev_buffers.cell_end,dev_buffers.density,dev_buffers.pressure,
                              cuda_constants);
        force_computation.Stop();

        cuMeasureNeighborsPerParticle(blocks,threads,N,
                                        dev_buffers.dev_rx,dev_buffers.dev_ry,
                                        dev_buffers.particle_cell_id,
                                        num_cells,nx,ny,
                                        dev_buffers.cell_start,dev_buffers.cell_end,
                                        cuda_constants);


        std::cout<<StepHeader(i,container,domain);


        cu_copy_buffer.Start();
        cuCopyBufferToHost(particles.fx.data(),dev_buffers.dev_fx,N);
        cuCopyBufferToHost(particles.fy.data(),dev_buffers.dev_fy,N);
        cuCopyBufferToHost(particles.rx.data(),dev_buffers.dev_rx,N);
        cuCopyBufferToHost(particles.ry.data(),dev_buffers.dev_ry,N);
        cuCopyBufferToHost(particles.vx.data(),dev_buffers.dev_vx,N);
        cuCopyBufferToHost(particles.vy.data(),dev_buffers.dev_vy,N);
        cu_copy_buffer.Stop();
    }   
    main_loop.Stop();

    std::cout<<"Main loop ="<<main_loop.GetTotalTimeInSeconds()<<" s"<<std::endl;
    std::cout<<force_computation.GetName()<<"="<<force_computation.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<density_approximation.GetName()<<"="<<density_approximation.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<container_build.GetName()<<"="<<container_build.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<cu_copy_buffer.GetName()<<"="<<cu_copy_buffer.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<integrate_operations.GetName()<<"="<<integrate_operations.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<checkpointing.GetName()<<"="<<checkpointing.GetTotalTimeInSeconds()<<"s"<<std::endl;
    

    cuFreeMemory(dev_buffers.dev_rx);
    cuFreeMemory(dev_buffers.dev_ry);
    cuFreeMemory(dev_buffers.dev_fx);
    cuFreeMemory(dev_buffers.dev_fy);
    cuFreeMemory(dev_buffers.dev_vx);
    cuFreeMemory(dev_buffers.dev_vy);
    cuFreeMemory(dev_buffers.particle_cell_id);

    cuFreeMemory(dev_buffers.cell_start);
    cuFreeMemory(dev_buffers.cell_end);

    cuFreeMemory(dev_buffers.density);
    cuFreeMemory(dev_buffers.pressure);

    return 0;
}


std::string StepHeader(int step, Container& particles, Domain& domain){
    std::ostringstream out;

    out<<"Step= "<<step<<"\t N="<<ParticlesInDomain(domain,particles)<<"\t E_k="<<ComputeKineticEnergy(particles)<<"\t |v|_max="<<ComputeMaxVelocity(particles)<<std::endl;

    return out.str();
}