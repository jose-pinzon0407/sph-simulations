#pragma once

#include "Container.h"
#include "SoAParticle.h"

class SoAContainer:public Container{

    protected:

    Particles particles;
    size_t size;

    public: 
    Particles& GetParticles(){
        return particles;
    }
    Particle operator[](const std::size_t i){
        return Particle{particles.density[i],particles.mass[i],pos(i),vel(i),force(i)};
    }

    void Build(){

    }

    void Rebuild(){

    }

    Vec3d pos(size_t i);
    Vec3d vel(size_t i);
    Vec3d force(size_t i);

    uint GetParticleIndex(size_t i);
    void SetParticleIndex(size_t i, uint idx);
    uint GetParticleCellIndex(size_t i);
    void SetParticleCellIndex(size_t i, uint idx);


    Real GetMass(size_t i);
    void SetDensity(size_t i, Real rho_i);
    Real GetDensity(size_t i);
    ScalarField& GetDensity(){

    }

    void SetForce(size_t i, Real fix, Real fiy, Real fiz);
    void SetForce(size_t i, Vec3d F);
    void AddForce(std::size_t i, Vec3d F);
    void AddForce(std::size_t i, Real fix, Real fiy, Real fiz);

    void SetPosition(size_t i, size_t j, Real val);
    void SetPosition(size_t, Real fix, Real fiy, Real fiz);
    void SetPosition(size_t, Vec3d R);
    void UpdatePosition(std::size_t i, Real dt);

    void SetVelocity(size_t i, size_t j, Real val);
    void SetVelocity(size_t, Real fix, Real fiy, Real fiz);
    void SetVelocity(size_t, Vec3d R);
    void UpdateVelocity(std::size_t i, Real dt);


    void Resize(size_t N) override;
    size_t Size();
    


};