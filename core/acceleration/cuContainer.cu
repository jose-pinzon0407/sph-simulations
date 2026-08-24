#include "cudaLaunchers.cuh"



void cuDeviceProperties(){
    int number_devices{0};
    cudaGetDeviceCount(&number_devices);
    std::cout<<"Total devices="<<number_devices<<std::endl;
    for(int i=0;i<number_devices;++i){
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        std::cout<<"[Context] Device name: "<<prop.name<<std::endl;
        std::cout<<"[Context] Max. threads per block: "<<prop.maxThreadsPerBlock<<std::endl;
        std::cout<<"[Context] Max. threads per multi.proc.: "<<prop.maxThreadsPerMultiProcessor<<std::endl;
        std::cout<<"[Context] Max. threads per dimension per block: "<<prop.maxThreadsDim[0]<<","<<prop.maxThreadsDim[1]<<","<<prop.maxThreadsDim[2]<<std::endl;
        std::cout<<"[Context] Max. blocks per multiprocessor: "<<prop.maxBlocksPerMultiProcessor<<std::endl;
        std::cout<<"[Context] Max. grid size: "<<prop.maxGridSize[0]<<","<<prop.maxGridSize[1]<<","<<prop.maxGridSize[2]<<std::endl;
        std::cout<<"[Context] Max. global memory [bytes]: "<<(double)prop.totalGlobalMem/1000000000.0<<"[GB]"<<std::endl;
        std::cout<<"[Context] In doubles: "<<(double)prop.totalGlobalMem/sizeof(double)/1.0e6<<" million doubles"<<std::endl;
        std::cout<<"[Context] Shared memory per SM: "<<(double)prop.sharedMemPerMultiprocessor/1e6<<"[MB]"<<std::endl;
        std::cout<<"[Context] Shared memory per block: "<<(double)prop.sharedMemPerBlock/1e6<<"[MB]"<<std::endl;
        std::cout<<"[Context] Registers per SM: "<<prop.regsPerMultiprocessor<<std::endl;
        std::cout<<"[Context] Registers per block: "<<prop.regsPerBlock<<std::endl;
        std::cout<<"[Context] Total const memory [bytes]: "<<prop.totalConstMem<<std::endl;
        std::cout<<"[Context] L2 cache size[bytes]: "<<prop.l2CacheSize<<std::endl;
    }
    
    cudaDeviceSynchronize();
}

void cuCreateContainer(SoAContainer& container, cuParticles& dev_particles){
    int N = container.Size();

    cuSetupBuffer(N,dev_particles.dev_rx);
    cuSetupBuffer(N,dev_particles.dev_ry);

    cuSetupBuffer(N, dev_particles.dev_vx);
    cuSetupBuffer(N, dev_particles.dev_vy);

    cuSetupBuffer(N, dev_particles.dev_fx);
    cuSetupBuffer(N, dev_particles.dev_fy);

    cuSetupBuffer(N, dev_particles.density);
    cuSetupBuffer(N, dev_particles.pressure);

    cuSetupBuffer<int>(N, dev_particles.particle_cell_id);

}

void cuReadContainer(SoAContainer& container, cuParticles& dev_particles){

    auto& particles = container.GetParticles();
    int N = container.Size();

    //copy positions
    cuCopyBufferToDevice(N,dev_particles.dev_rx,particles.rx.data());
    cuCopyBufferToDevice(N,dev_particles.dev_ry,particles.ry.data());

    //copy velocities
    cuCopyBufferToDevice(N,dev_particles.dev_vx,particles.vx.data());
    cuCopyBufferToDevice(N,dev_particles.dev_vy,particles.vy.data());

    // cuSetBufferValue<Real>(N,dev_particles.dev_vx,0);
    // cuSetBufferValue<Real>(N,dev_particles.dev_vy,0);

     cudaDeviceSynchronize();

}

void cuWriteContainer(SoAContainer& container, cuParticles& dev_particles){
    int N = container.Size();
    auto& particles = container.GetParticles();
    // cudaDeviceSynchronize();
    //copy positions
    cuCopyBufferToHost(N,particles.rx.data(),dev_particles.dev_rx);
    cuCopyBufferToHost(N,particles.ry.data(),dev_particles.dev_ry);

    //copy velocities
    cuCopyBufferToHost(N,particles.vx.data(),dev_particles.dev_vx);
    cuCopyBufferToHost(N,particles.vy.data(),dev_particles.dev_vy);

    //copy forces
    cuCopyBufferToHost(N,particles.fx.data(),dev_particles.dev_fx);
    cuCopyBufferToHost(N,particles.fy.data(),dev_particles.dev_fy);
    cudaDeviceSynchronize();
    

}

void cuSetupNeighborLists(SoAContainer& container, cuNeighborList& neighbor_list){
    int N = container.Size();

    cuSetupBuffer((N+1),neighbor_list.dev_offset_per_particle);
    cuSetupBuffer(N, neighbor_list.dev_neighbors_per_particle);
    cudaDeviceSynchronize();

    std::cout<<"[cuConfiguration]Configured buffers for neighbor lists\n";

}