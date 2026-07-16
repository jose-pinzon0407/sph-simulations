#include "container/SortedContainer.h"

void SortedSoAContainer::Build(){
    std::sort(particles.particle_id.begin(),particles.particle_id.end(),[&](uint a, uint b){
        return particles.cell_id[a]<particles.cell_id[b];
    }
    );

    SortIntoCells();
}

void SortedSoAContainer::SortIntoCells(){
    ScalarField ordered_density;
    ordered_density.resize(Size());



    std::vector<uint> cell_id_ordered;
    cell_id_ordered.resize(Size());

    std::vector<uint> particle_id_ordered;
    particle_id_ordered.resize(Size());
        
    std::vector<Real> ordered_rx;
    std::vector<Real> ordered_ry;
    std::vector<Real> ordered_rz;

    ordered_rx.resize(Size());
    ordered_ry.resize(Size());
    ordered_rz.resize(Size());


    std::vector<Real> ordered_vx;
    std::vector<Real> ordered_vy;
    std::vector<Real> ordered_vz;

    ordered_vx.resize(Size());
    ordered_vy.resize(Size());
    ordered_vz.resize(Size());

    std::vector<Real> ordered_fx;
    std::vector<Real> ordered_fy;
    std::vector<Real> ordered_fz;


    ordered_fx.resize(Size());
    ordered_fy.resize(Size());
    ordered_fz.resize(Size());

    for(size_t i=0;i<Size();++i){
            
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

void SortedSoAContainer::Rebuild(){
    
}

