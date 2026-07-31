#include "ListsForceComputation.h"

void ListsForceComputation::ApproximateDensity(SoAContainer& particles){
    Particles& soa_particles = particles.GetParticles();
    Real mu=MU;
    size_t N = particles.Size();

    for(size_t i=0;i<N;++i){
        std::vector<int>& i_neighbors = neighbor_list.GetNeighborsOfParticle(i);

        Real accum_rho_i=0.0;

        Real rx_i = soa_particles.rx[i];
        Real ry_i = soa_particles.ry[i];
        Real rz_i = soa_particles.rz[i];


        size_t num_neighbors = i_neighbors.size();
        for(int n=0;n<num_neighbors;++n){
            
            int j_idx = i_neighbors[n];

            Real dx = rx_i-soa_particles.rx[j_idx];
            Real dy = ry_i-soa_particles.ry[j_idx];
            Real dz = rz_i-soa_particles.rz[j_idx];

            Real r2 = dx*dx+dy*dy+dz*dz;

            accum_rho_i += DensityW(r2)*particles.GetMass(j_idx);

        }//end neighbor list iterator
        
        particles.SetDensity(i,accum_rho_i);
    }
}

void ListsForceComputation::CalculatePressure(SoAContainer& particles){
    size_t N = particles.Size();

    if(pressure_field.size() != N){
        pressure_field.resize(N,0.0);
    }

    for(size_t i=0;i<N;++i){
        pressure_field[i] = GAS_CONSTANT*(particles.GetDensity(i)-RHO_0);
    }
}

void ListsForceComputation::ComputeForces(SoAContainer& particles){
    
    ApproximateDensity(particles);
    CalculatePressure(particles);

    Particles& soa_particles = particles.GetParticles();
    Real mu=MU;
    size_t N = particles.Size();
    Real fg_y=-10.0;

    #pragma omp parallel for schedule(runtime)
    for(size_t i=0;i<N;++i){
        std::vector<int>& i_neighbors = neighbor_list.GetNeighborsOfParticle(i);

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


        size_t num_neighbors = i_neighbors.size();
        for(int n=0;n<num_neighbors;++n){
            
            int j_idx = i_neighbors[n];

            Real dx = rx_i-soa_particles.rx[j_idx];
            Real dy = ry_i-soa_particles.ry[j_idx];
            Real dz = rz_i-soa_particles.rz[j_idx];

            Real r2 = dx*dx+dy*dy+dz*dz;


            Real r1 = std::sqrt(r2);
            Real rho_j = soa_particles.density[j_idx];

            //Pressure forces
            if(r1> TOL){
            Real pressure_multiplier = -0.5*(pressure_field[i]/(rho_i*rho_i)+pressure_field[j_idx]/(rho_j*rho_j))*soa_particles.mass[j_idx];
            
            Real h_r12=(h1-r1)*(h1-r1);
            fp_x += -3.0*pressure_constant*h_r12/r1*dx*pressure_multiplier;
            fp_y += -3.0*pressure_constant*h_r12/r1*dy*pressure_multiplier;
            fp_z += -3.0*pressure_constant*h_r12/r1*dz*pressure_multiplier;
            
            }


            //Viscosity forces
            Real dvx = soa_particles.vx[j_idx]-vi_x;
            Real dvy = soa_particles.vy[j_idx]-vi_y;
            Real dvz = soa_particles.vz[j_idx]-vi_z;
            Real viscosity_multiplier = mu/(soa_particles.density[j_idx]*M_PI*h6)*45.0*(H-r1);

            fv_x += viscosity_multiplier*dvx;
            fv_y += viscosity_multiplier*dvy;
            fv_z += viscosity_multiplier*dvz;


        }//end neighbor list iterator
        
        soa_particles.fx[i]+=fp_x;
        soa_particles.fy[i]+=fp_y;
        soa_particles.fz[i]+=fp_z;

        soa_particles.fx[i]+=fv_x;
        soa_particles.fy[i]+=fv_y;
        soa_particles.fz[i]+=fv_z;
    
    }
}