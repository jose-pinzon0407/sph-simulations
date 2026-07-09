#include "SoAContainer.h"

size_t SoAContainer::Size(){
    return size;
}

void SoAContainer::Resize(size_t N){
    this->size = N;

    particles.mass.resize(N,1.0);
    particles.density.resize(N);
    particles.fx.resize(N);
    particles.fy.resize(N);
    particles.fz.resize(N);

    particles.rx.resize(N);
    particles.ry.resize(N);
    particles.rz.resize(N);

    particles.vx.resize(N);
    particles.vy.resize(N);
    particles.vz.resize(N);

    particles.cell_id.resize(N);
    particles.particle_id.resize(N);

}

Vec3d SoAContainer::pos(size_t i){
    return Vec3d{particles.rx[i],particles.ry[i],particles.rz[i]};
}

Vec3d SoAContainer::vel(size_t i){
    return Vec3d{particles.vx[i],particles.vy[i],particles.vz[i]};
}

Vec3d SoAContainer::force(size_t i){
    return {particles.fx[i],particles.fy[i],particles.fz[i]};
}

uint SoAContainer::GetParticleIndex(size_t i){
    return particles.particle_id[i];
}

void SoAContainer::SetParticleIndex(size_t i, uint idx){
    particles.particle_id[i]=idx;
}

uint SoAContainer::GetParticleCellIndex(size_t i){
    return particles.cell_id[i];
}

void SoAContainer::SetParticleCellIndex(size_t i, uint idx){
    particles.cell_id[i]=idx;
}

Real SoAContainer::GetMass(size_t i){
    return particles.mass[i];
}
void SoAContainer::SetDensity(size_t i, Real rho_i){
    particles.density[i]=rho_i;
}
Real SoAContainer::GetDensity(size_t i){
    return particles.density[i];
}


void SoAContainer::SetForce(size_t i, Real fix, Real fiy, Real fiz){
    particles.fx[i]=fix;
    particles.fy[i]=fiy;
    particles.fz[i]=fiz;
}
void SoAContainer::SetForce(size_t i, Vec3d F){
    SetForce(i,F[0],F[1],F[2]);
}
void SoAContainer::AddForce(std::size_t i, Vec3d F){
    AddForce(i,F[0],F[1],F[2]);
}
void SoAContainer::AddForce(std::size_t i, Real fix, Real fiy, Real fiz){
    particles.fx[i]+=fix;
    particles.fy[i]+=fiy;
    particles.fz[i]+=fiz;
}

void SoAContainer::SetPosition(size_t i, size_t j, Real val){
    if(j==0){
        particles.rx[i]=val;
        return;
    }

    if(j==1){
        particles.ry[i]=val;
        return;
    }

    if(j==2){
        particles.rz[i]=val;
        return;
    }
}
void SoAContainer::SetPosition(size_t i, Real rix, Real riy, Real riz){
    particles.rx[i]=rix;
    particles.ry[i]=riy;
    particles.rz[i]=riz;
}
void SoAContainer::SetPosition(size_t i, Vec3d R){
    SetPosition(i, R[0],R[1],R[2]);
}
void SoAContainer::UpdatePosition(size_t i, Real dt){
    particles.rx[i] += dt*particles.vx[i];
    particles.ry[i] += dt*particles.vy[i];
    particles.rz[i] += dt*particles.vz[i];
}

void SoAContainer::SetVelocity(size_t i, size_t j, Real val){
    if(j==0) particles.vx[i]=val;
    if(j==1) particles.vy[i]=val;
    if(j==2) particles.vz[i]=val;
}
void SoAContainer::SetVelocity(size_t i, Real vix, Real viy, Real viz){
    particles.vx[i]=vix;
    particles.vy[i]=viy;
    particles.vz[i]=viz;
}
void SoAContainer::SetVelocity(size_t i, Vec3d V){
    SetVelocity(i, V[0],V[1],V[2]);
}
void SoAContainer::UpdateVelocity(std::size_t i, Real dt){
    particles.vx[i] += dt*particles.fx[i];
    particles.vy[i] += dt*particles.fy[i];
    particles.vz[i] += dt*particles.fz[i];
}