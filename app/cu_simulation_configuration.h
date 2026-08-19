#pragma once

#include"acceleration/cudaLaunchers.cuh"

/**
 * Simulation launchers for CUDA enabled simulation
*/

#include "simulation_configuration.h"
using namespace simulation;
namespace gpu_simulation{

    cuGridInformation grid;
    cuParticles dev_particles;
    cuNeighborList dev_neighbor_lists;
    int blocks, threads;
    Real* inv_rho;


    void cudaConfigure(Domain& domain, Real cs, int thrds, int N){
        grid.cell_size=cs;
        Vec3d domain_lenghts=domain.GetLength();
        grid.nx= std::ceil(domain_lenghts[0]/grid.cell_size);
        grid.ny= std::ceil(domain_lenghts[1]/grid.cell_size);
        grid.nz= std::max((int)std::ceil(domain_lenghts[2]/grid.cell_size),1);


        cuSetupBuffer(grid.GetTotalCells(),dev_particles.cell_start);
        cuSetupBuffer(grid.GetTotalCells(),dev_particles.cell_end);

        std::cout<<"[cuInitializeGrid] Using uniform grid with: "<<grid.GetTotalCells()<<" cells."<<std::endl;
        std::cout<<"[cuInitializeGrid] Cells per direction: "<<grid.nx<<","<<grid.ny<<","<<grid.nz<<std::endl;
        std::cout<<"[cuInitializeGrid] Cell size: "<<grid.cell_size<<std::endl;

        threads=thrds;
        blocks = (N+threads-1)/threads;
        std::cout<<"[GPU] Launch configuration: "<<blocks<<" blocks with "<<threads<<" threads grid\n";
        std::cout<<"[GPU] Total threads: "<<blocks*threads<<std::endl;
    }

    void cudaSimulationSetup(SoAContainer& container){
        cuDeviceProperties();
        cuCreateContainer(container,dev_particles);

        cuSetupBuffer(container.Size(),inv_rho);

        cuReadContainer(container,dev_particles);
        cuCreateCellBuffers(grid.GetTotalCells(),dev_particles);
        cuSetupNeighborLists(container,dev_neighbor_lists);
    }


    void PreLoopComputationWithCellNeighbors(SoAContainer& container,cuConstants& cuda_constants){
        int N=container.Size();


        timers.Start("neighbor_search");
        cuKernelFindNeighbors(blocks,threads,N,
                              dev_particles.particle_cell_id,
                              grid.nx,grid.ny,grid.cell_size,
                              dev_particles.cell_start,dev_particles.cell_end,
                              dev_particles.dev_rx,dev_particles.dev_ry,
                              dev_particles.dev_vx,dev_particles.dev_vy,
                              dev_particles.dev_fx,dev_particles.dev_fy
                            );
        timers.Stop("neighbor_search");

        std::vector<int> host_cell_id;
        std::vector<int> host_cell_start;
        std::vector<int> host_cell_end;
        host_cell_id.resize(N);
        host_cell_start.resize(grid.GetTotalCells());
        host_cell_end.resize(grid.GetTotalCells());
        cuCopyBufferToHost(N,host_cell_id.data(),dev_particles.particle_cell_id);
        cuCopyBufferToHost(grid.GetTotalCells(),host_cell_start.data(),dev_particles.cell_start);
        cuCopyBufferToHost(grid.GetTotalCells(),host_cell_end.data(),  dev_particles.cell_end);
     
        PrintCellsToFile(0,host_cell_start,host_cell_end);

        timers.Start("approximate_density");
        cuKernelApproximateDensitySorted(blocks,threads,N,dev_particles.dev_rx,dev_particles.dev_ry,dev_particles.particle_cell_id,grid.GetTotalCells(),grid.nx,grid.ny,dev_particles.cell_start,dev_particles.cell_end,dev_particles.density,cuda_constants);
        timers.Stop("approximate_density");

        timers.Start("buffer_inversion");
        cuInvertBuffer(N,blocks,threads,dev_particles.density,inv_rho);
        timers.Stop("buffer_inversion");


        timers.Start("calculate_pressure");
        cuKernelComputePressureSorted(blocks,threads,N,dev_particles.density,dev_particles.pressure);
        timers.Stop("calculate_pressure");


        std::vector<Real> host_pressure;
        host_pressure.resize(N);

        timers.Start("force_computation");
        cuKernelComputeForcesSorted(blocks,threads,N,
            dev_particles.dev_rx,dev_particles.dev_ry,
            dev_particles.dev_fx,dev_particles.dev_fy,
            dev_particles.dev_vx,dev_particles.dev_vy,
            dev_particles.particle_cell_id,
            grid.GetTotalCells(),grid.nx,grid.ny,
            dev_particles.cell_start,dev_particles.cell_end,
            inv_rho, dev_particles.pressure,
            cuda_constants);
        timers.Stop("force_computation");

        timers.Start("cuda_buffer_copy");
        cuWriteContainer(container,dev_particles);
        timers.Stop("cuda_buffer_copy");


        std::cout<<"[Simulation] Output step 0 positions\n";
        checkpoint->WriteCheckPoint(0,container);

        std::cout<<StepHeader(0,container,domain);


    }

    void MaiLoopComputationWithCellNeighbors(SoAContainer& container,cuConstants& cuda_constants){
        int N=container.Size();
        Integrator integrator{dt,&domain};

        timers.Start("main_loop");
        for(int i=1;i<steps;++i){
            
            integrator.UpdateVelocities(container);
            integrator.UpdatePositions(container);
            checkpoint->WriteCheckPoint(i,container);

            timers.Start("cuda_buffer_copy");
            cuReadContainer(container,dev_particles);
            timers.Stop("cuda_buffer_copy");

            timers.Start("neighbor_search");
            cuKernelFindNeighbors(blocks,threads,
                                 container.Size(),
                                 dev_particles.particle_cell_id,
                                 grid.nx,grid.ny,grid.cell_size,
                                 dev_particles.cell_start,dev_particles.cell_end,
                                 dev_particles.dev_rx,dev_particles.dev_ry,
                                 dev_particles.dev_vx,dev_particles.dev_vy,
                                 dev_particles.dev_fx,dev_particles.dev_fy
                                );
            timers.Stop("neighbor_search");


            std::vector<int> host_cell_id;
            std::vector<int> host_cell_start;
            std::vector<int> host_cell_end;
            host_cell_id.resize(N);
            host_cell_start.resize(grid.GetTotalCells());
            host_cell_end.resize(grid.GetTotalCells());
            cuCopyBufferToHost(N,host_cell_id.data(),dev_particles.particle_cell_id);
            cuCopyBufferToHost(grid.GetTotalCells(),host_cell_start.data(),dev_particles.cell_start);
            cuCopyBufferToHost(grid.GetTotalCells(),host_cell_end.data(),  dev_particles.cell_end);
     
            PrintCellsToFile(i,host_cell_start,host_cell_end);



            timers.Start("approximate_density");
            cuKernelApproximateDensitySorted(blocks,threads,N,
                dev_particles.dev_rx,dev_particles.dev_ry,
                dev_particles.particle_cell_id,
                grid.GetTotalCells(),grid.nx,grid.ny,
                dev_particles.cell_start,dev_particles.cell_end,
                dev_particles.density,cuda_constants);
            timers.Stop("approximate_density");
            
            timers.Start("buffer_inversion");
            cuInvertBuffer(N,blocks,threads,dev_particles.density,inv_rho);
            timers.Stop("buffer_inversion");

            timers.Start("calculate_pressure");
            cuKernelComputePressureSorted(blocks,threads,N,dev_particles.density,dev_particles.pressure);
            timers.Stop("calculate_pressure");

            timers.Start("force_computation");
            cuKernelComputeForcesSorted(blocks,threads,N,
                dev_particles.dev_rx,dev_particles.dev_ry,
                dev_particles.dev_fx,dev_particles.dev_fy,
                dev_particles.dev_vx,dev_particles.dev_vy,
                dev_particles.particle_cell_id,
                grid.GetTotalCells(),grid.nx,grid.ny,
                dev_particles.cell_start,dev_particles.cell_end,
                inv_rho,dev_particles.pressure,
                cuda_constants);
            timers.Stop("force_computation");


            std::cout<<StepHeader(i,container,domain);

            timers.Start("cuda_buffer_copy");
            cuWriteContainer(container,dev_particles);
            timers.Stop("cuda_buffer_copy");



        }
        timers.Stop("main_loop");


    }

    void RunSimulationWithCellNeighbors(){
        cuConstants cuda_constants{H,MU};

        PreLoopComputationWithCellNeighbors(simulation::container,cuda_constants);
        MaiLoopComputationWithCellNeighbors(simulation::container,cuda_constants);

    }


    /*****************************************************************************
     * LISTS LAUNCHERS
     ****************************************************************************/

    void PreLoopComputationWithLists(SoAContainer& container,cuConstants& cuda_constants){
        int N=container.Size();


        timers.Start("neighbor_search");
        cuKernelFindNeighbors(blocks,threads,N,
                              dev_particles.particle_cell_id,
                              grid.nx,grid.ny,grid.cell_size,
                              dev_particles.cell_start,dev_particles.cell_end,
                              dev_particles.dev_rx,dev_particles.dev_ry,
                              dev_particles.dev_vx,dev_particles.dev_vy,
                              dev_particles.dev_fx,dev_particles.dev_fy
                            );
        timers.Stop("neighbor_search");



        cudaMemset(dev_neighbor_lists.dev_neighbors_per_particle,0,N*sizeof(int));
        cudaMemset(dev_neighbor_lists.dev_offset_per_particle,0,(N+1)*sizeof(int));

        timers.Start("neighbor_list");
        cuKernelNeighborLists(blocks,threads,
                                  N,
                                  dev_particles.dev_rx,dev_particles.dev_ry,
                                  dev_particles.particle_cell_id,
                                  grid.GetTotalCells(),grid.nx,grid.ny,
                                  dev_particles.cell_start,dev_particles.cell_end,
                                  dev_neighbor_lists.dev_neighbors_per_particle,dev_neighbor_lists.dev_offset_per_particle, dev_neighbor_lists.dev_neighbor_indices
                                  );
        timers.Stop("neighbor_list");


        timers.Start("approximate_density");
        cuKernelApproximateDensityWithNeighborLists(blocks,threads,
                                                    N,
                                                    dev_particles.dev_rx,dev_particles.dev_ry,
                                                    dev_particles.particle_cell_id,
                                                    grid.GetTotalCells(),grid.nx,grid.ny,
                                                    dev_particles.density,
                                                    dev_neighbor_lists.dev_offset_per_particle,dev_neighbor_lists.dev_neighbor_indices,
                                                    cuda_constants);
        timers.Stop("approximate_density");

        timers.Start("buffer_inversion");
        cuInvertBuffer(N,blocks,threads,dev_particles.density,inv_rho);
        timers.Stop("buffer_inversion");

        timers.Start("calculate_pressure");
        cuKernelComputePressureSorted(blocks,threads,N,dev_particles.density,dev_particles.pressure);
        timers.Stop("calculate_pressure");


        timers.Start("force_computation");
        cuKernelComputeForcesWithNeighborLists(blocks,threads,N,
            dev_particles.dev_rx,dev_particles.dev_ry,
            dev_particles.dev_fx,dev_particles.dev_fy,
            dev_particles.dev_vx,dev_particles.dev_vy,
            dev_particles.particle_cell_id,
            grid.GetTotalCells(),grid.nx,grid.ny,
            inv_rho,dev_particles.pressure,
            dev_neighbor_lists.dev_offset_per_particle,dev_neighbor_lists.dev_neighbor_indices,
            cuda_constants);
        timers.Stop("force_computation");

        timers.Start("cuda_buffer_copy");
        cuWriteContainer(container,dev_particles);
        timers.Stop("cuda_buffer_copy");


        std::cout<<"[Simulation] Output step 0 positions\n";
        checkpoint->WriteCheckPoint(0,container);

        std::cout<<StepHeader(0,container,domain);


    }

    void MaiLoopComputationWithLists(SoAContainer& container,cuConstants& cuda_constants){
        int N=container.Size();
        Integrator integrator{dt,&domain};

        timers.Start("main_loop");
        for(int i=1;i<steps;++i){
            
            integrator.UpdateVelocities(container);
            integrator.UpdatePositions(container);
            checkpoint->WriteCheckPoint(i,container);

            timers.Start("cuda_buffer_copy");
            cuReadContainer(container,dev_particles);
            timers.Stop("cuda_buffer_copy");

            timers.Start("neighbor_search");
            cuKernelFindNeighbors(blocks,threads,
                                 container.Size(),
                                 dev_particles.particle_cell_id,
                                 grid.nx,grid.ny,grid.cell_size,
                                 dev_particles.cell_start,dev_particles.cell_end,
                                 dev_particles.dev_rx,dev_particles.dev_ry,
                                 dev_particles.dev_vx,dev_particles.dev_vy,
                                 dev_particles.dev_fx,dev_particles.dev_fy
                                );
            timers.Stop("neighbor_search");
            
            cudaMemset(dev_neighbor_lists.dev_neighbors_per_particle,0,N*sizeof(int));
            cudaMemset(dev_neighbor_lists.dev_offset_per_particle,0,(N+1)*sizeof(int));

            timers.Start("neighbor_list");
            cuKernelNeighborLists(blocks,threads,
                                  N,
                                  dev_particles.dev_rx,dev_particles.dev_ry,
                                  dev_particles.particle_cell_id,
                                  grid.GetTotalCells(),grid.nx,grid.ny,
                                  dev_particles.cell_start,dev_particles.cell_end,
                                  dev_neighbor_lists.dev_neighbors_per_particle,dev_neighbor_lists.dev_offset_per_particle, dev_neighbor_lists.dev_neighbor_indices
                                  );
            timers.Stop("neighbor_list");

        std::vector<int> host_neighbors_per_particle;
        std::vector<int> host_offset_per_particle;
        std::vector<int> host_indices_per_particle;
        
        host_neighbors_per_particle.resize(N);
           host_offset_per_particle.resize(N+1);
          host_indices_per_particle.resize(N);

        cuCopyBufferToHost(N,host_indices_per_particle.data(),dev_neighbor_lists.dev_neighbor_indices);
        cuCopyBufferToHost(N,host_neighbors_per_particle.data(),dev_neighbor_lists.dev_neighbors_per_particle);
        cuCopyBufferToHost(N+1,host_offset_per_particle.data(),dev_neighbor_lists.dev_offset_per_particle);

            timers.Start("approximate_density");
            cuKernelApproximateDensityWithNeighborLists(blocks,threads,N,
                dev_particles.dev_rx,dev_particles.dev_ry,
                dev_particles.particle_cell_id,
                grid.GetTotalCells(),grid.nx,grid.ny,
                dev_particles.density,
                dev_neighbor_lists.dev_offset_per_particle,dev_neighbor_lists.dev_neighbor_indices,
                cuda_constants);
            timers.Stop("approximate_density");

            timers.Start("buffer_inversion");
            cuInvertBuffer(N,blocks,threads,dev_particles.density,inv_rho);
            timers.Stop("buffer_inversion");



            timers.Start("calculate_pressure");
            cuKernelComputePressureSorted(blocks,threads,N,dev_particles.density,dev_particles.pressure);
            timers.Stop("calculate_pressure");



            timers.Start("force_computation");
            cuKernelComputeForcesWithNeighborLists(blocks,threads,N,
                dev_particles.dev_rx,dev_particles.dev_ry,
                dev_particles.dev_fx,dev_particles.dev_fy,
                dev_particles.dev_vx,dev_particles.dev_vy,
                dev_particles.particle_cell_id,
                grid.GetTotalCells(),grid.nx,grid.ny,
                inv_rho,dev_particles.pressure,
                dev_neighbor_lists.dev_offset_per_particle,dev_neighbor_lists.dev_neighbor_indices,
                cuda_constants);
            timers.Stop("force_computation");

            // timers.Start("force_computation");
            // cuKernelComputeForcesSorted(blocks,threads,N,
            //     dev_particles.dev_rx,dev_particles.dev_ry,
            //     dev_particles.dev_fx,dev_particles.dev_fy,
            //     dev_particles.dev_vx,dev_particles.dev_vy,
            //     dev_particles.particle_cell_id,
            //     grid.GetTotalCells(),grid.nx,grid.ny,
            //     dev_particles.cell_start,dev_particles.cell_end,
            //     dev_particles.density,dev_particles.pressure,
            //     cuda_constants);
            // timers.Stop("force_computation");
            

            std::cout<<StepHeader(i,container,domain);

            timers.Start("cuda_buffer_copy");
            cuWriteContainer(container,dev_particles);
            timers.Stop("cuda_buffer_copy");



        }
        timers.Stop("main_loop");


    }

    void RunSimulationWithLists(){
        cuConstants cuda_constants{H,MU};


        PreLoopComputationWithLists(simulation::container,cuda_constants);

        MaiLoopComputationWithLists(simulation::container,cuda_constants);

    }



    void RunSimulationGPU(){
        std::cout<<"[Configuration] Launching CUDA-accelerated simulation\n";
        
        if(parser.GetParameter("compute_type")=="mp-lists"){
            std::cout<<"[Configuration] Running simulation with lists\n";
            RunSimulationWithLists();
            return;
        }

        if(parser.GetParameter("compute_type")=="sorted"){
            std::cout<<"[Configuration] Running simulation with sorted cells\n";
            RunSimulationWithCellNeighbors();
            return;
        }

    }

}


// #endif