#pragma once

#include "Particle.h"

class Container{
    public:

    virtual Particle operator[](const std::size_t i)=0;


    virtual void Build()=0;
    virtual void Rebuild()=0;

    virtual Vec3d pos(size_t i)=0;
    virtual Vec3d vel(size_t i)=0;
    virtual Vec3d force(size_t i)=0;

    virtual size_t Size()=0;
    virtual uint GetParticleIndex(size_t i)=0;
    virtual void SetParticleIndex(size_t i, uint idx)=0;
    virtual uint GetParticleCellIndex(size_t i)=0;
    virtual void SetParticleCellIndex(size_t i, uint idx)=0;

    virtual void Resize(std::size_t N)=0;
    virtual Real GetMass(size_t i)=0;
    virtual void SetDensity(size_t i, Real rho)=0;
    virtual Real GetDensity(size_t)=0;
    virtual ScalarField& GetDensity()=0;
    virtual void SetForce(size_t i, Real fix, Real fiy, Real fiz)=0;
    virtual void SetForce(size_t i, Vec3d F)=0;
    virtual void AddForce(std::size_t i, Vec3d F)=0;
    virtual void AddForce(std::size_t i, Real fix, Real fiy, Real fiz)=0;

    virtual void SetPosition(size_t i, size_t j, Real val)=0;
    virtual void SetPosition(size_t, Real fix, Real fiy, Real fiz)=0;
    virtual void SetPosition(size_t, Vec3d R)=0;
    virtual void UpdatePosition(std::size_t i, Real dt)=0;
    virtual void SetVelocity(size_t i, size_t j, Real val)=0;
    virtual void SetVelocity(size_t, Real fix, Real fiy, Real fiz)=0;
    virtual void SetVelocity(size_t, Vec3d R)=0;
    virtual void UpdateVelocity(std::size_t i, Real dt)=0;
};

class ParticleContainer:public Container{

    public:
    ParticleContainer(size_t N);
    Particle operator[](const std::size_t i) override;

    void Build()override{

    }
    void Rebuild() override{

    }

    Vec3d pos(size_t i)override;
    Vec3d vel(size_t i)override;
    Vec3d force(size_t i)override;

    size_t Size()override;
    uint GetParticleIndex(size_t i)override;
    void SetParticleIndex(size_t i, uint idx)override;
    uint GetParticleCellIndex(size_t i)override;
    void SetParticleCellIndex(size_t i, uint idx)override;

    
    Real GetMass(size_t i)override;
    void SetDensity(size_t i, Real rho)override;
    Real GetDensity(size_t)override;
    ScalarField& GetDensity()override;
    void SetForce(size_t i, Real fix, Real fiy, Real fiz)override;
    void SetForce(size_t i, Vec3d F)override;
    void AddForce(std::size_t i, Vec3d F)override;
    void AddForce(std::size_t i, Real fix, Real fiy, Real fiz)override;

    void SetPosition(size_t i, size_t j, Real val)override;
    void SetPosition(size_t, Real fix, Real fiy, Real fiz)override;
    void SetPosition(size_t, Vec3d R)override;
    void UpdatePosition(std::size_t i, Real dt)override;
    void SetVelocity(size_t i, size_t j, Real val)override;
    void SetVelocity(size_t, Real fix, Real fiy, Real fiz)override;
    void SetVelocity(size_t, Vec3d R)override;
    void UpdateVelocity(std::size_t i, Real dt)override;
    void Resize(std::size_t N)override;

    protected:
    std::vector<Particle> particles;

};