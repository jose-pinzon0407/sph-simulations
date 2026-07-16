#include "SoAForceComputation.h"


void SoAForceComputation::ComputeForces(SoAContainer& particles){

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

void SoAForceComputation::ApproximateDensity(SoAContainer& particles){
    size_t N = particles.Size();

    for(size_t i=0;i<N;++i){
        Real accum_rho_i=0.0;
        for(size_t j=0;j<N;++j){
            Real dx = particles.GetParticles().rx[i]-particles.GetParticles().rx[j];
            Real dy = particles.GetParticles().ry[i]-particles.GetParticles().ry[j];
            Real dz = particles.GetParticles().rz[i]-particles.GetParticles().rz[j];

            Real r2 = dx*dx+dy*dy+dz*dz;

            if(r2 > h2)continue;

            accum_rho_i += DensityW(r2)*particles.GetMass(j);

        }
        particles.SetDensity(i,accum_rho_i);
    }
}

void SoAForceComputation::CalculatePressure(SoAContainer& particles){

    size_t N = particles.Size();

    if(pressure_field.size() != N){
        pressure_field.resize(N,0.0);
    }

    for(size_t i=0;i<particles.Size();++i){
        pressure_field[i] = GAS_CONSTANT*(particles.GetDensity(i)-RHO_0);
    }
}