#include "Particle.h"

Particle::Particle(Real dens, Real ms, Vec3d& pos, Vec3d& vel, Vec3d& force):rho{dens},mass{ms},r{pos},v{vel},f{force}{

}

Real Particle::GetRho(){
    return rho;
}

Real Particle::GetMass(){
    return mass;
}

Vec3d& Particle::GetR(){
    return r;
}

Vec3d& Particle::GetV(){
    return v;
}

Vec3d& Particle::GetF(){
    return f;
}