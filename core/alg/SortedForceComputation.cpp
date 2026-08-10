#include "SoAForceComputation.h"

void SortedForceComputation::ApproximateDensity(SoAContainer& particles){

    size_t N = particles.Size();
    // std::vector<int> sorted_neighbor_counts;
    // sorted_neighbor_counts.resize(N);
    // std::fill(sorted_neighbor_counts.begin(),sorted_neighbor_counts.end(),0);
    for(size_t i=0;i<N;++i){
        Real accum_rho_i=0.0;

        uint cell_i = particles.GetParticleCellIndex(i);
    
        // int cell_start = neighbor_search.GetCellStart(cell_i);
        // int cell_end = neighbor_search.GetCellEnd(cell_i);


        //neighboring cells
        Vec3i multi_index_i = neighbor_search.Unflatten(cell_i);



        if(neighbor_search.Flatten(multi_index_i) != cell_i){
            std::cout<<"Error in flattening\n";
            std::exit(EXIT_FAILURE);
        }

        if(neighbor_search.Flatten(neighbor_search.Unflatten(cell_i))!=particles.GetParticleCellIndex(i)){
            std::cout<<"flatten(unflatten) issue\n";
            std::exit(EXIT_FAILURE);
        }

        for(int c1=-1;c1<=1;++c1){
            for(int c2=-1;c2<=1;++c2){

                Vec3i multi_index_j = multi_index_i;// + Vec3i{c1,c2,0};

                multi_index_j[0] += c1;

                if(multi_index_j[0] < 0 || multi_index_j[0] >= neighbor_search.GetCellsInDirectionX()) continue;

                multi_index_j[1] += c2;
                if(multi_index_j[1] < 0|| multi_index_j[1] >= neighbor_search.GetCellsInDirectionY()) continue;


                // if(multi_index_j[0] < 0|| multi_index_j[0] >= neighbor_search.GetCellsInDirectionX()) continue;

                int cell_j = neighbor_search.Flatten(multi_index_j);


                if(neighbor_search.Unflatten(cell_j) != multi_index_j){
                    std::cout<<"Error in flattening\n";
                    std::exit(EXIT_FAILURE);
                }

                if(cell_j >= neighbor_search.GetCellsInDirectionX()*neighbor_search.GetCellsInDirectionY()){
                    std::cout<<"Invalid global index\n";
                    std::exit(EXIT_FAILURE);
                }

                int cell_start_j = neighbor_search.GetCellStart(cell_j);
                int cell_end_j = neighbor_search.GetCellEnd(cell_j);

                if(cell_start_j < 0 ||  cell_end_j<0){
                    continue;
                }

                for(int pj=cell_start_j;pj<cell_end_j;++pj){

                     if(neighbor_search.Flatten(neighbor_search.Unflatten(cell_j))!=particles.GetParticleCellIndex(pj)){
                    std::cout<<"flatten(unflatten) issue\n";
                    std::exit(EXIT_FAILURE);
                    }

                    Real dx = particles.GetParticles().rx[i]-particles.GetParticles().rx[pj];
                    Real dy = particles.GetParticles().ry[i]-particles.GetParticles().ry[pj];
                    Real dz = particles.GetParticles().rz[i]-particles.GetParticles().rz[pj];

                    Real r2 = dx*dx+dy*dy+dz*dz;


                    if(r2 > h2)continue;

                    // sorted_neighbor_counts[i]++;
                    accum_rho_i += DensityW(r2)*particles.GetMass(pj);

                }//end iterate pj cell_j


            }
        }

        particles.SetDensity(i,accum_rho_i);

    }//end i iteration
}


void SortedForceComputation::CalculatePressure(SoAContainer& particles){

    size_t N = particles.Size();

    if(pressure_field.size() != N){
        pressure_field.resize(N,0.0);
    }

    for(size_t i=0;i<particles.Size();++i){
        pressure_field[i] = GAS_CONSTANT*(particles.GetDensity(i)-RHO_0);
    }
}

void SortedForceComputation::ComputeForces(SoAContainer& particles){
    Real mu=MU;

    ApproximateDensity(particles);
    CalculatePressure(particles);

    Real fg_y=-10.0;
    size_t N = particles.Size();

    Particles& soa_particles = particles.GetParticles();

    #pragma omp parallel for schedule(runtime)
    for(size_t i=0;i<N;++i){

        uint cell_i = particles.GetParticleCellIndex(i);
        Vec3i multi_index_i = neighbor_search.Unflatten(cell_i);


        //Pressure force
        Real fp_x=0.0;
        Real fp_y=0.0;
        Real fp_z=0.0;

        //Viscosity force
        Real fv_x=0.0;
        Real fv_y=0.0;
        Real fv_z=0.0;

        soa_particles.fx[i]=0.0;
        soa_particles.fy[i]=fg_y;
        soa_particles.fz[i]=0.0;

        Real vi_x = soa_particles.vx[i];
        Real vi_y = soa_particles.vy[i];
        Real vi_z = soa_particles.vz[i];

        Real rho_i = soa_particles.density[i];

        Real rx_i = soa_particles.rx[i];
        Real ry_i = soa_particles.ry[i];
        Real rz_i = soa_particles.rz[i];

        for(int cx=-1;cx<=1;++cx){
            for(int cy=-1;cy<=1;++cy){

                Vec3i multi_index_j = multi_index_i;
                
                //iterate neighbors
                multi_index_j[0]+=cx;

                if(multi_index_j[0] < 0 || multi_index_j[0] >= neighbor_search.GetCellsInDirectionX()) continue;

                multi_index_j[1]+=cy;
                if(multi_index_j[1] < 0 || multi_index_j[1] >= neighbor_search.GetCellsInDirectionY()) continue;

                int cell_j = neighbor_search.Flatten(multi_index_j);

                int cell_start_j = neighbor_search.GetCellStart(cell_j);
                int cell_end_j=neighbor_search.GetCellEnd(cell_j);

                if(cell_start_j<0 || cell_end_j<0){
                    continue;
                }

                for(int pj=cell_start_j;pj<cell_end_j;++pj){

                    Real dx = soa_particles.rx[i]-soa_particles.rx[pj];
                    Real dy = soa_particles.ry[i]-soa_particles.ry[pj];
                    Real dz = soa_particles.rz[i]-soa_particles.rz[pj];

                    Real r2 = dx*dx+dy*dy+dz*dz;

                    if(r2 > h2)continue;

                    Real r1 = std::sqrt(r2);
                    Real rho_j = soa_particles.density[pj];

                    //Pressure forces
                    if(r1> TOL){
                        Real pressure_multiplier = -0.5*(pressure_field[i]/(rho_i*rho_i)+pressure_field[pj]/(rho_j*rho_j))*soa_particles.mass[pj];
                // Real pressure_multiplier = -0.5*(pressure_field[i]/(rho_i*rho_i)+pressure_field[j]/(rho_j*rho_j));

                        Real h_r12=(h1-r1)*(h1-r1);
                        fp_x += -3.0*pressure_constant*h_r12/r1*dx*pressure_multiplier;
                        fp_y += -3.0*pressure_constant*h_r12/r1*dy*pressure_multiplier;
                        fp_z += -3.0*pressure_constant*h_r12/r1*dz*pressure_multiplier;
            
            
                    }

                    //Viscosity forces
                    Real dvx = soa_particles.vx[pj]-vi_x;
                    Real dvy = soa_particles.vy[pj]-vi_y;
                    Real dvz = soa_particles.vz[pj]-vi_z;
                    Real viscosity_multiplier = mu/(soa_particles.density[pj]*M_PI*h6)*45.0*(H-r1);

                    fv_x += viscosity_multiplier*dvx;
                    fv_y += viscosity_multiplier*dvy;
                    fv_z += viscosity_multiplier*dvz;
                }//end cell_j particle iterator



            }
        }

        soa_particles.fx[i]+=fp_x;
        soa_particles.fy[i]+=fp_y;
        soa_particles.fz[i]+=fp_z;

        soa_particles.fx[i]+=fv_x;
        soa_particles.fy[i]+=fv_y;
        soa_particles.fz[i]+=fv_z;

    }//end i-iterate

}
