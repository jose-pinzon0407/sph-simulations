#include "Integrator.h"

Integrator::Integrator(Real ss,Domain* dom):dt{ss},domain{dom}{

}

void Integrator::UpdatePositions(Container& particles){
    size_t N = particles.Size();

    for(size_t i=0;i<N;++i){
        particles.UpdatePosition(i,dt);

        Vec3d r_i = particles.pos(i);
        Vec3d v_i = particles.vel(i);
        ApplyBoundaryConditions(r_i,v_i);
        particles.SetPosition(i,r_i);
        particles.SetVelocity(i,v_i);

    }

}

void Integrator::UpdateVelocities(Container& particles){
    size_t N = particles.Size();
    for(size_t i=0;i<N;++i){
        particles.UpdateVelocity(i,dt);
    }
}

bool Integrator::ApplyBoundaryConditions(Vec3d& r, Vec3d& v){
    bool apply=false;
        for(size_t i=0;i<r.size();++i){
            if(r[i]>domain->h_corner[i]){
                apply=true;
                r[i] = domain->h_corner[i];
                v[i] *= -1;
            }

            if(r[i]<domain->l_corner[i]){
                apply=true;
                r[i] = domain->l_corner[i];
                v[i] *=-1;   
            }
        }

        return apply;
}