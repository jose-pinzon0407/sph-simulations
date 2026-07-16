#pragma once

#include "container/Container.h"
#include "Domain.h"
Real ComputeKineticEnergy(Container& particles){
    Real e_k=0.0;
    for(size_t i=0;i<particles.Size();++i){
        e_k += 0.5*SquaredNorm(particles[i].v)*particles[i].mass;
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