#include "SoAForceComputation.h"

void SortedForceComputation::ApproximateDensity(SortedSoAContainer& particles){

    size_t N = particles.Size();
    std::vector<int> sorted_neighbor_counts;
    sorted_neighbor_counts.resize(N);
    std::fill(sorted_neighbor_counts.begin(),sorted_neighbor_counts.end(),0);
    for(size_t i=0;i<N;++i){
        Real accum_rho_i=0.0;

        uint cell_i = particles.GetParticleCellIndex(i);
    
        int cell_start = neighbor_search.GetCellStart(cell_i);
        int cell_end = neighbor_search.GetCellEnd(cell_i);

        // // self-cell
        // for(int j = cell_start;j<cell_end;++j){
        //     Real dx = particles.GetParticles().rx[i]-particles.GetParticles().rx[j];
        //     Real dy = particles.GetParticles().ry[i]-particles.GetParticles().ry[j];
        //     Real dz = particles.GetParticles().rz[i]-particles.GetParticles().rz[j];

        //     Real r2 = dx*dx+dy*dy+dz*dz;

        //     if(r2 > h2)continue;

        //     accum_rho_i += DensityW(r2)*particles.GetMass(j);

        // }


        //neighboring cells
        Vec3i multi_index_i = neighbor_search.Unflatten(cell_i);

        //   if(cell_i == 100 && i == 1359){
        //     std::cout<<"Cell 100 particle 1359\n";
        //     std::cout<<"i: "<<i<<std::endl;
        //     std::cout<<"cell i: "<<cell_i<<std::endl;
        //     std::cout<<"local i: "<<multi_index_i<<std::endl;
        //     std::cout<<"i_rx:"<<particles.GetParticles().rx[i]<<std::endl;
        //     std::cout<<"i_ry:"<<particles.GetParticles().ry[i]<<std::endl;
        //     std::cout<<"i_rz:"<<particles.GetParticles().rz[i]<<std::endl;

        // }

        if(neighbor_search.Flatten(multi_index_i) != cell_i){
            std::cout<<"Error in flattening\n";
            std::exit(EXIT_FAILURE);
        }

        if(neighbor_search.Flatten(neighbor_search.Unflatten(cell_i))!=particles.GetParticleCellIndex(i)){
            std::cout<<"flatten(unflatten) issue\n";
            std::exit(EXIT_FAILURE);
        }

        // std::cout<<"Particle: "<<i<<" in cell: "<<multi_index_i<<" rx:"<<particles.GetParticles().rx[i]<<" ry:"<<particles.GetParticles().ry[i]<<std::endl;
        // return;
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

                // std::cout<<"Cell i: "<<multi_index_i<<" interacts with j:"<<multi_index_j<<std::endl;


                // if(cell_j == 100 && cell_i==100){
                //     std::cout<<"Begin:"<<cell_start_j<<std::endl;
                //     std::cout<<"End:"<<cell_end_j<<std::endl;
                // }
                for(int pj=cell_start_j;pj<cell_end_j;++pj){

                     if(neighbor_search.Flatten(neighbor_search.Unflatten(cell_j))!=particles.GetParticleCellIndex(pj)){
                    std::cout<<"flatten(unflatten) issue\n";
                    std::exit(EXIT_FAILURE);
                    }

                    Real dx = particles.GetParticles().rx[i]-particles.GetParticles().rx[pj];
                    Real dy = particles.GetParticles().ry[i]-particles.GetParticles().ry[pj];
                    Real dz = particles.GetParticles().rz[i]-particles.GetParticles().rz[pj];

                    Real r2 = dx*dx+dy*dy+dz*dz;

                    // if((cell_j==100 && cell_i==100) && i == 1359 ){
                    //     std::cout<<"Cell i: "<<cell_i<<" cell j: "<<cell_j<<" i particle : "<<i<<std::endl;
                    //     std::cout<<"particle j:"<<pj<<std::endl;
                    //     std::cout<<"local j: "<<multi_index_j<<std::endl;
                    //     std::cout<<"local i: "<<multi_index_i<<std::endl;
                    //     std::cout<<"pj_rx:"<<particles.GetParticles().rx[pj]<<std::endl;
                    //     std::cout<<"pj_ry:"<<particles.GetParticles().ry[pj]<<std::endl;
                    //     std::cout<<"pj_rz:"<<particles.GetParticles().rz[pj]<<std::endl;

                    //     std::cout<<"r2:"<<dx*dx+dy*dy+dz*dz<<std::endl;

                    // }

                    if(r2 > h2)continue;
                    // std::cout<<dx<<" "<<dy<<" "<<dz<<"\n";

                    sorted_neighbor_counts[i]++;
                    accum_rho_i += DensityW(r2)*particles.GetMass(pj);

                }//end iterate pj cell_j


            }
        }

        particles.SetDensity(i,accum_rho_i);

    }//end i iteration


    for(size_t i=0;i<N;++i){
        int neighbors_i=0;
        Real accum_rho_i=0.0;
        for(size_t j=0;j<N;++j){
            Real dx = particles.GetParticles().rx[i]-particles.GetParticles().rx[j];
            Real dy = particles.GetParticles().ry[i]-particles.GetParticles().ry[j];
            Real dz = particles.GetParticles().rz[i]-particles.GetParticles().rz[j];

            Real r2 = dx*dx+dy*dy+dz*dz;

            if(r2 > h2)continue;

            accum_rho_i += DensityW(r2)*particles.GetMass(j);
            neighbors_i++;
        }
        Real current_rho_i = particles.GetDensity(i);
        if( std::abs(current_rho_i - accum_rho_i) > 1e-7 ){
            std::cout<<"Denstiy divergence\n";
            std::cout<<"rho_sorted:"<<current_rho_i<<"- accum_rho: "<<accum_rho_i<<std::endl;
            std::cout<<"sorted neighbors:"<<sorted_neighbor_counts[i]<<"- nosort neighbors: "<<neighbors_i<<std::endl;
            std::cout<<"i:"<<i<<" rx:"<<particles.GetParticles().rx[i]<<" ry: "<<particles.GetParticles().ry[i]<<" cell:"<<particles.GetParticleCellIndex(i)<<std::endl;
            std::exit(EXIT_FAILURE);
        }
        // particles.SetDensity(i,accum_rho_i);
    }
}


void SortedForceComputation::CalculatePressure(SortedSoAContainer& particles){

    size_t N = particles.Size();

    if(pressure_field.size() != N){
        pressure_field.resize(N,0.0);
    }

    for(size_t i=0;i<particles.Size();++i){
        pressure_field[i] = GAS_CONSTANT*(particles.GetDensity(i)-RHO_0);
    }
}

void SortedForceComputation::ComputeForces(SortedSoAContainer& particles){

    Real mu = MU;

    ApproximateDensity(particles);
    CalculatePressure(particles);
    Real fg_y = -10.0;
    size_t N = particles.Size();

    Particles& soa_particles = particles.GetParticles();

    for(size_t i=0;i<N;++i){

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

        for(size_t j=0;j<N;++j){

            Real dx = soa_particles.rx[i]-soa_particles.rx[j];
            Real dy = soa_particles.ry[i]-soa_particles.ry[j];
            Real dz = soa_particles.rz[i]-soa_particles.rz[j];

            Real r2 = dx*dx+dy*dy+dz*dz;

            if(r2 > h2)continue;

            Real r1 = std::sqrt(r2);

            Real rho_j = soa_particles.density[j];

            //Forces of pressure kernel
            if(r1>TOL){
                Real pressure_multiplier = -0.5*(pressure_field[i]/(rho_i*rho_i)+pressure_field[j]/(rho_j*rho_j))*soa_particles.mass[j];
                // Real pressure_multiplier = -0.5*(pressure_field[i]/(rho_i*rho_i)+pressure_field[j]/(rho_j*rho_j));

                Real h_r12=(h1-r1)*(h1-r1);
                fp_x += -3.0*pressure_constant*h_r12/r1*dx*pressure_multiplier;
                fp_y += -3.0*pressure_constant*h_r12/r1*dy*pressure_multiplier;
                fp_z += -3.0*pressure_constant*h_r12/r1*dz*pressure_multiplier;
            
            
            }
            


            //Forces of viscosity kernel

            Real dvx = soa_particles.vx[j]-vi_x;
            Real dvy = soa_particles.vy[j]-vi_y;
            Real dvz = soa_particles.vz[j]-vi_z;
            Real viscosity_multiplier = mu/(soa_particles.density[j]*M_PI*h6)*45.0*(H-r1);

            fv_x += viscosity_multiplier*dvx;
            fv_y += viscosity_multiplier*dvy;
            fv_z += viscosity_multiplier*dvz;



        }


        soa_particles.fx[i]+=fp_x;
        soa_particles.fy[i]+=fp_y;
        soa_particles.fz[i]+=fp_z;

        soa_particles.fx[i]+=fv_x;
        soa_particles.fy[i]+=fv_y;
        soa_particles.fz[i]+=fv_z;

    }

}