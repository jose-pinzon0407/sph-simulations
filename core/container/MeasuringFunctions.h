#pragma once

#include "container/Container.h"
#include "Domain.h"
Real ComputeKineticEnergy(Container& particles){
    Real e_k=0.0;
    for(size_t i=0;i<particles.Size();++i){

        Vec3d v_i = particles[i].v;
        if(!std::isfinite(v_i[0])
           || !std::isfinite(v_i[1])
           || !std::isfinite(v_i[2])
            ){
                std::cout<<"Invalid velocity at particle "<<i<<": "<<v_i<<"\n";
                break;            
            }
        if(!std::isfinite(particles[i].mass)
            ){
                std::cout<<"Invalid mass at particle "<<i<<": "<<"\n";
                break;            
            }

        e_k += 0.5*SquaredNorm(particles[i].v)*particles[i].mass;

        if(!std::isfinite(e_k)
            ){
                std::cout<<"Invalid EK  "<<i<<": "<<e_k<<"\n";
                break;            
            }
            
            Real my_ek = 0.5*SquaredNorm(particles[i].v)*particles[i].mass;

        if(!std::isfinite(my_ek)
            ){
                std::cout<<"Invalid EK at particle "<<i<<": "<<my_ek<<"\n";
                break;            
            }
    }

    return e_k;
}

Real ComputeMaxVelocity(Container& particles){
    Real v_max=0.0;

    for(size_t i=0;i<particles.Size();++i){

        v_max = std::max((Real)SquaredNorm(particles[i].v),v_max);
    }
    return v_max;
}

int ParticlesInDomain(Domain& domain, Container& particles){
    int total_particles=0;
    for(size_t i=0;i<particles.Size();++i){
        if(domain.IsInside(particles[i].r)){
            total_particles++;
        }
    }

    return total_particles;
}