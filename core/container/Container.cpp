#include "Container.h"


// Particle ParticleContainer::operator[](const size_t i){
//     return particles[i];
// }

void ParticleContainer::Resize(size_t N){
    particles.resize(N);
}

// size_t ParticleContainer::Size(){
//     return particles.size();
// }

// Vec3d ParticleContainer::pos(size_t i){
//     return particles[i].r;
// }

// Vec3d ParticleContainer::vel(size_t i){
//     return particles[i].v;
// }

Vec3d ParticleContainer::force(size_t i) {
    return particles[i].f;
}

uint ParticleContainer::GetParticleIndex(size_t i){
    return i;
}

void ParticleContainer::SetParticleIndex(size_t i, uint idx){
    
}

uint ParticleContainer::GetParticleCellIndex(size_t i){
    return i;
}

void ParticleContainer::SetParticleCellIndex(size_t i, uint idx){

}

// Real ParticleContainer::GetMass(size_t i){
//     return particles[i].mass;
// }

// Real ParticleContainer::GetDensity(size_t i){
//     return particles[i].rho;
// }

ScalarField& ParticleContainer::GetDensity(){
 //not implemented
 
}

// void ParticleContainer::SetDensity(size_t i, Real dens){
//     particles[i].rho=dens;
// }

// void ParticleContainer::SetForce(size_t i, Real fix, Real fiy, Real fiz){
//     particles[i].f[0]=fix;
//     particles[i].f[1]=fiy;
//     particles[i].f[2]=fiz;
// }

// void ParticleContainer::SetForce(size_t i, Vec3d F){
//     SetForce(i,F[0],F[1],F[2]);
// }

// void ParticleContainer::AddForce(size_t i, Vec3d F){
//     AddForce(i,F[0],F[1],F[2]);

// }

// void ParticleContainer::AddForce(size_t i, Real fix, Real fiy, Real fiz){
//     particles[i].f[0]+=fix;
//     particles[i].f[1]+=fiy;
//     particles[i].f[2]+=fiz;
// }


void ParticleContainer::SetPosition(size_t i, Real rix, Real riy, Real riz){
    particles[i].r[0]=rix;
    particles[i].r[1]=riy;
    particles[i].r[2]=riz;
}

void ParticleContainer::SetPosition(size_t i, Vec3d F){
    SetPosition(i,F[0],F[1],F[2]);
}

void ParticleContainer::SetPosition(size_t i, size_t j, Real val){
    particles[i].r[j]=val;
}

void ParticleContainer::UpdatePosition(size_t i, Real dt){
    particles[i].r[0] += dt*particles[i].v[0];
    particles[i].r[1] += dt*particles[i].v[1];
    particles[i].r[2] += dt*particles[i].v[2];
}

void ParticleContainer::SetVelocity(size_t i, Real vix, Real viy, Real viz){
    particles[i].v[0]=vix;
    particles[i].v[1]=viy;
    particles[i].v[2]=viz;
}

void ParticleContainer::SetVelocity(size_t i, Vec3d V){
    SetVelocity(i,V[0],V[1],V[2]);
}

void ParticleContainer::SetVelocity(size_t i, size_t j, Real val){
    particles[i].v[j]=val;
}


void ParticleContainer::UpdateVelocity(size_t i, Real dt){
    particles[i].v[0] += dt*particles[i].f[0];
    particles[i].v[1] += dt*particles[i].f[1];
    particles[i].v[2] += dt*particles[i].f[2];
}