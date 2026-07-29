#include<thrust/sort.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>

#include<numeric>
#include<algorithm>
#include "CUDAKernels.cuh"

__global__
void cuSortParticlesIntoCells(int N, Real* rx, Real* ry, Real cell_size, int* particle_cell_index, int nx, int ny);

void cuThrustCallBack(int N, int* particle_cell_index, Real* rx, Real* ry 
                                                     ,Real* fx, Real* fy
                                                     ,Real* vx, Real* vy);

__global__
void cuSetCellRanges(int N, int* cell_start, int* cell_end, int* particle_cell_index);

__global__
void cuNeighborsPerParticle(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, cuConstants kernel_constants, int* neighs, int* cneighs);

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

void cuMeasureNeighborsPerParticle(int blocks, int threads, int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, cuConstants kernel_constants){

    thrust::host_vector<int> h_neighbors_pp(N,0);
    thrust::host_vector<int> h_cneighbors_pp(N,0);
    thrust::device_vector<int> neighbors_per_particle=h_neighbors_pp;
    thrust::device_vector<int> candidate_neighbors_per_particle=h_neighbors_pp;

    cuNeighborsPerParticle<<<blocks,threads>>>(N,rx,ry,particle_cell_id,num_cells,nx,ny,start,end,kernel_constants,thrust::raw_pointer_cast(neighbors_per_particle.data()),thrust::raw_pointer_cast(candidate_neighbors_per_particle.data()));
    cudaDeviceSynchronize();

    h_neighbors_pp = neighbors_per_particle;
    h_cneighbors_pp = candidate_neighbors_per_particle;
    int average = std::reduce(h_neighbors_pp.begin(),h_neighbors_pp.end());
    auto mx = std::max_element(h_neighbors_pp.begin(),h_neighbors_pp.end());
    auto mn = std::min_element(h_neighbors_pp.begin(),h_neighbors_pp.end());

    double avg = average/N;



    std::cout<<"[HOST] Average neighbors per particle= "<<avg<<std::endl;
    std::cout<<"[HOST] Max. element= "<<*mx<<std::endl;
    std::cout<<"[HOST] Min. element= "<<*mn<<std::endl;

    double std_dev=0.0;
    for(int i=0;i<N;++i){
        std_dev += (h_neighbors_pp[i]-avg)*(h_neighbors_pp[i]-avg);
    }  
    
    std::cout<<"[HOST] Std. dev= "<<sqrt(std_dev/N)<<std::endl;

    double caverage = std::reduce(h_cneighbors_pp.begin(),h_cneighbors_pp.end())/N;

    std::cout<<"[HOST] Averga candidate neighbors per particle= "<<caverage<<std::endl;
    mx = std::max_element(h_cneighbors_pp.begin(),h_cneighbors_pp.end());
    mn = std::min_element(h_cneighbors_pp.begin(),h_cneighbors_pp.end());
    std::cout<<"[HOST] Max. candidate element= "<<*mx<<std::endl;
    std::cout<<"[HOST] Min. candidate element= "<<*mn<<std::endl;

    std_dev=0.0;
    for(int i=0;i<N;++i){
        std_dev += (h_cneighbors_pp[i]-caverage)*(h_cneighbors_pp[i]-caverage);
    }  
    
    std::cout<<"[HOST] Std. candidate dev= "<<sqrt(std_dev/N)<<std::endl;

}

/***
 * 
 * KERNEL DEFINITIONS
 * 
 * 
 */

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

void cuOutputDeviceProperties(){
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


__global__
void cuNeighborsPerParticle(int N, Real* rx, Real* ry, int* particle_cell_id, int num_cells, int nx, int ny, int* start, int* end, cuConstants kernel_constants, int* neighbors, int* cneighs){

    int idx=threadIdx.x+blockIdx.x*blockDim.x;

    if(idx>=N) return;


    Real rix = rx[idx];
    Real riy = ry[idx];


 


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

                cneighs[idx]++;
                              
                if(r2>kernel_constants.h2)continue;

                neighbors[idx]++;



            }

        }
    }


}