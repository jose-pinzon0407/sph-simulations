#include "ForceComputation.h"

ForceComputation::ForceComputation(Real r,size_t num_parts):kernels{r},h2{r*r}{
    pressure_field.resize(num_parts);
    constant = 315.0/(64.0*M_PI*std::pow(r,9.0));
}


void ForceComputation::ComputeForces(ParticleContainer& particles){
    Real mu=MU;

    ApproximateDensity(particles);
    CalculatePressure(particles);
    
    Vec3d force_gravity=Vec3d{0,-10,0};

    size_t N = particles.Size();
    for(size_t i=0;i<N;++i){
        Vec3d r_i = particles.pos(i);
        Vec3d force_pressure_i{0,0,0};
        Vec3d force_viscosity_i{0,0,0};

        Real rho_i = particles.GetDensity(i);
        Real pressure_i = pressure_field[i];
        particles.SetForce(i,force_gravity);

        for(size_t j=0;j<N;++j){
            // Vec3d rij = particles[i].r-particles[j].r;
            Vec3d rij = r_i-particles.pos(j);
            if(SquaredNorm(rij)>h2) continue;

            Real rho_j = particles.GetDensity(j);
            Real pressure_j = pressure_field[j];
            Real pressure_multiplier = -0.5*(pressure_i/(rho_i*rho_i)+pressure_j/(rho_j*rho_j))*particles[j].mass;
            Vec3d pressure_gradient = kernels.pressure.GradW(rij);

            force_pressure_i += pressure_multiplier*pressure_gradient;

            Real viscosity_multiplier = mu/rho_j*kernels.viscosity.LaplacianW(rij);
            Vec3d velocity_difference = particles.vel(j)-particles.vel(i);

            force_viscosity_i += viscosity_multiplier*velocity_difference;

        }

        particles.AddForce(i,force_pressure_i);
        particles.AddForce(i,force_viscosity_i);
    }


}

void ForceComputation::CalculatePressure(ParticleContainer& particles){
    std::fill(pressure_field.begin(),pressure_field.end(),0.0);
    for(int i=0;i<particles.Size();++i){
        pressure_field[i] = GAS_CONSTANT*(particles.GetDensity(i)-RHO_0);
    }
}

void ForceComputation::ApproximateDensity(ParticleContainer& particles){
    const size_t n = particles.Size();
    for(int i=0;i<n;++i){
        Real accum_rho_i=0.0;
        Vec3d r_i = particles.pos(i);
        for(int j=0;j<n;++j){

            // Vec3d rij = particles[i].r-particles[j].r;
            Vec3d rij = r_i-particles[j].r;
            // Vec3d rij = r_i-particles.pos(j);
            double r2 = SquaredNorm(rij);
            if(r2>h2) continue;

            Real mass_j = particles.GetMass(j);
            accum_rho_i += kernels.density.W(rij)*mass_j;
            // accum_rho_i += (h2-r2)*(h2-r2)*(h2-r2)*constant * mass_j;
        }

        particles.SetDensity(i,accum_rho_i);
    }
}