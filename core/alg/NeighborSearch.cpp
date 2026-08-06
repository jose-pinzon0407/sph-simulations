#include "NeighborSearch.h"


NeighborSearch::NeighborSearch(Domain& domain, Real cell_size){
    SetGrid(domain,cell_size);
}

void NeighborSearch::FindNeighbors(SoAContainer& particles){
    size_t num_particles = particles.Size();

    SetParticleCellIds(particles);
    Sort(particles);
    SetCellRanges(num_particles,particles);
    // OutputOrder(particles);
    
}


void NeighborSearch::SetParticleCellIds(SoAContainer& container){
    size_t num_particles=container.Size();
    auto& particles = container.GetParticles();

    #pragma omp parallel for schedule(runtime)
    for(size_t i=0;i<num_particles;++i){
        Real rxi = particles.rx[i];
        Real ryi = particles.ry[i];
        Real rzi = particles.rz[i];

        if(rxi < 0.0) std::cout<<"Out of domain\n";
        if(ryi < 0.0) std::cout<<"Out of domain\n";
        if(rzi < 0.0) std::cout<<"Out of domain\n";

        int ci0 = std::floor(rxi/grid.cell_size);
        int ci1 = std::floor(ryi/grid.cell_size);
        int ci2 = std::floor(rzi/grid.cell_size);

        ci0 = std::min(ci0, grid.nx-1);
        ci1 = std::min(ci1, grid.ny-1);
        ci2 = std::min(ci2,grid.nz-1);

        // int ci0 = (rxi/grid.cell_size);
        // int ci1 = (ryi/grid.cell_size);
        // int ci2 = (rzi/grid.cell_size);

        uint idx = ci0 + ci1*grid.nx;
        particles.cell_id[i]=idx;
    }


}

void NeighborSearch::OutputOrder(SoAContainer& particles){
    std::ofstream output;            
    std::string name = "order_step_"+std::to_string(step)+".txt";
    output.open(name.c_str());
    for(size_t i=0;i<particles.Size();++i){
        output<<i<<"\t"<<particles.GetParticleIndex(i)<<"\t"<<particles.GetParticleCellIndex(i)<<"\t";
        output<<particles.pos(i)[0]<<"\t"<<particles.pos(i)[1];
        output<<std::endl;
    }

    output.close();
    name = "cell_order_step_"+std::to_string(step)+".txt";
    output.open(name.c_str());

    for(size_t i=0;i<grid.GetTotalCells();++i){
        output<<"Cell "<<i<<"["<<grid.cell_start[i]<<","<<grid.cell_end[i]<<"]"<<std::endl;
    }

    output.close();

    step++;
}

void NeighborSearch::SetGrid(Domain& domain, Real radius){
        Vec3d box_length = domain.GetLength();

        grid.cell_size = radius;
        grid.nx = std::ceil(box_length[0]/grid.cell_size);
        grid.ny = std::ceil(box_length[1]/grid.cell_size);
        grid.nz = std::max((int)std::ceil(box_length[2]/grid.cell_size),1);

        grid.ResizeCellRanges();

        std::cout<<"[NeighborSearch] Using uniform grid with: "<<grid.GetTotalCells()<<" cells."<<std::endl;
        std::cout<<"[NeighborSearch] Cells per direction: "<<grid.nx<<","<<grid.ny<<","<<grid.nz<<std::endl;
        std::cout<<"[NeighborSearch] Cell size: "<<grid.cell_size<<std::endl;
}

void NeighborSearch::SetCellRanges(uint N, SoAContainer& particles){
    grid.RefillCellRanges();
    for(uint i=0;i<N;++i){
        uint cell_i_id = particles.GetParticleCellIndex(i);
        if(i==0 || cell_i_id != particles.GetParticleCellIndex(i-1)){
            grid.cell_start[cell_i_id] = i;
        }

        if(i==(N-1) || cell_i_id != particles.GetParticleCellIndex(i+1)){
            grid.cell_end[cell_i_id]= i+1;//open interval
        }
    }


    SortValidation(N);
}

void NeighborSearch::SortValidation(uint N){
    int accum_particles=0;
    for(size_t i=0;i<grid.GetTotalCells();++i){
        if(grid.cell_end[i] > -1 && grid.cell_start[i]> -1)
            accum_particles += (grid.cell_end[i]-grid.cell_start[i]);
    }
        // std::cout<<"[NeighborSearch] Validation, counted particles "<<accum_particles<<" vs. N="<<N<<std::endl;
}

void NeighborSearch::SortParticleIds(SoAContainer& container){
    auto& particles = container.GetParticles();

    std::sort(particles.particle_id.begin(),particles.particle_id.end(),[&](uint a, uint b){
        return particles.cell_id[a]<particles.cell_id[b];
              }
    );

}

void NeighborSearch::SortContainer(SoAContainer& container){
    int N = container.Size();
    ScalarField ordered_density;
    ordered_density.resize(N);



    std::vector<uint> cell_id_ordered;
    cell_id_ordered.resize(N);

    std::vector<uint> particle_id_ordered;
    particle_id_ordered.resize(N);
        
    std::vector<Real> ordered_rx;
    std::vector<Real> ordered_ry;
    std::vector<Real> ordered_rz;

    ordered_rx.resize(N);
    ordered_ry.resize(N);
    ordered_rz.resize(N);


    std::vector<Real> ordered_vx;
    std::vector<Real> ordered_vy;
    std::vector<Real> ordered_vz;

    ordered_vx.resize(N);
    ordered_vy.resize(N);
    ordered_vz.resize(N);

    std::vector<Real> ordered_fx;
    std::vector<Real> ordered_fy;
    std::vector<Real> ordered_fz;


    ordered_fx.resize(N);
    ordered_fy.resize(N);
    ordered_fz.resize(N);

    auto& particles = container.GetParticles();
    
    for(size_t i=0;i<N;++i){
                
        uint old_particle_index = particles.particle_id[i];
        //reorder rho
        ordered_density[i] = particles.density[old_particle_index];
        //reorder r
        ordered_rx[i] = particles.rx[old_particle_index];
        ordered_ry[i] = particles.ry[old_particle_index];
        ordered_rz[i] = particles.rz[old_particle_index];
        //reorder v
        ordered_vx[i] = particles.vx[old_particle_index];
        ordered_vy[i] = particles.vy[old_particle_index];
        ordered_vz[i] = particles.vz[old_particle_index];
        //reorder f
        ordered_fx[i] = particles.fx[old_particle_index];
        ordered_fy[i] = particles.fy[old_particle_index];
        ordered_fz[i] = particles.fz[old_particle_index];
        //reorder cell_id
        cell_id_ordered[i] = particles.cell_id[old_particle_index];
        //reorder particle id
        particles.particle_id[i]=i;

    }

    //assign vectors
    particles.density = ordered_density;

    particles.rx = ordered_rx;
    particles.ry = ordered_ry;
    particles.rz = ordered_rz;

    particles.vx = ordered_vx;
    particles.vy = ordered_vy;
    particles.vz = ordered_vz;
            
    particles.fx = ordered_fx;
    particles.fy = ordered_fy;
    particles.fz = ordered_fz;
    particles.cell_id=cell_id_ordered;

}

void NeighborSearch::Sort(SoAContainer& container){
    SortParticleIds(container);
    SortContainer(container);
}


std::vector<int>& NeighborSearch::GetEnd(){
    return grid.cell_end;
}

std::vector<int>& NeighborSearch::GetStart(){
    return grid.cell_start;
}

int NeighborSearch::GetCellStart(int cell){
    return grid.cell_start[cell];
}

int NeighborSearch::GetCellEnd(int cell){
    return grid.cell_end[cell];
}

int NeighborSearch::Flatten(const Vec3i& c){
    return grid.MapThreeToOne(c);
}

Vec3i NeighborSearch::Unflatten(unsigned int index){
    return grid.MapOneToThree(index);
}

int NeighborSearch::GetCellsInDirectionX(){
    return grid.nx;
}

int NeighborSearch::GetCellsInDirectionY(){
    return grid.ny;
}

int NeighborSearch::GetCellsInDirectionZ(){
    return grid.nz;
}

int NeighborSearch::GetTotalCells(){
    return grid.GetTotalCells();
}
