#pragma once

#include "Particle.h"
#include "SoAParticle.h"


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

    inline Particle operator[](const std::size_t i) override{
        return particles[i];
    }

    void Build()override{

    }
    void Rebuild() override{

    }

    inline Vec3d pos(size_t i)override{
        return particles[i].r;
    }

    inline Vec3d vel(size_t i)override{
        return particles[i].v;
    }

    inline Vec3d force(size_t i)override;

    inline size_t Size()override{
        return particles.size();
    }
    uint GetParticleIndex(size_t i)override;
    void SetParticleIndex(size_t i, uint idx)override;
    uint GetParticleCellIndex(size_t i)override;
    void SetParticleCellIndex(size_t i, uint idx)override;

    
    inline Real GetMass(size_t i)override{
        return particles[i].mass;
    }
    inline void SetDensity(size_t i, Real dens)override{
        particles[i].rho=dens;
    }
    inline Real GetDensity(size_t i)override{
        return particles[i].rho;
    }

    inline ScalarField& GetDensity()override;
    inline void SetForce(size_t i, Real fix, Real fiy, Real fiz)override{
        particles[i].f[0]=fix;
        particles[i].f[1]=fiy;
        particles[i].f[2]=fiz;
    }
    inline void SetForce(size_t i, Vec3d F)override{
        SetForce(i,F[0],F[1],F[2]);
    }
    
    inline void AddForce(std::size_t i, Vec3d F)override{
        AddForce(i,F[0],F[1],F[2]);
    }
    inline void AddForce(std::size_t i, Real fix, Real fiy, Real fiz)override{
        particles[i].f[0]+=fix;
        particles[i].f[1]+=fiy;
        particles[i].f[2]+=fiz;
    }

    inline void SetPosition(size_t i, size_t j, Real val)override;
    inline void SetPosition(size_t, Real fix, Real fiy, Real fiz)override;
    inline void SetPosition(size_t, Vec3d R)override;
    inline void UpdatePosition(std::size_t i, Real dt)override;
    inline void SetVelocity(size_t i, size_t j, Real val)override;
    inline void SetVelocity(size_t, Real fix, Real fiy, Real fiz)override;
    inline void SetVelocity(size_t, Vec3d R)override;
    inline void UpdateVelocity(std::size_t i, Real dt)override;
    inline void Resize(std::size_t N)override;

    protected:
    std::vector<Particle> particles;

};