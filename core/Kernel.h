#pragma once

#include "math/Math.h"

/**
 * Does the
 */

class SmoothingFunction{

    public:
    SmoothingFunction(Real r):h{r},h2{r*r}{

    }

    virtual Real W(Vec3d& rij)=0;
    virtual Vec3d GradW(Vec3d& r)=0;
    virtual Real LaplacianW(Vec3d& r)=0;


    protected:
    virtual void SetVal()=0;


    protected:
    Real h;
    Real h2;
    Real val=0.0;
};

class DensityKernel:public SmoothingFunction{
    public:
    DensityKernel(Real r);

    inline Real W(Vec3d& rij) override{
        Real r2 = SquaredNorm(rij);
        return  (h2-r2)*(h2-r2)*(h2-r2)*val;
    }

    inline Vec3d GradW(Vec3d& r)override;
    inline Real LaplacianW(Vec3d& r) override;

    protected:
    void SetVal()override;
};

class PressureKernel:public SmoothingFunction{
    public:
    PressureKernel(Real rij);
    inline Real W(Vec3d& rij)override;
    inline Vec3d GradW(Vec3d& rij)override{
        Real r1 = VecNorm(rij);
        if(r1<TOL) return Vec3d{0.0,0.0,0.0};

        Real dW= -3.0*val*(h-r1)*(h-r1);
        return (dW/r1)*rij;
    }
    inline Real LaplacianW(Vec3d& r) override;

    protected:
    void SetVal()override;
};

class ViscosityKernel:public SmoothingFunction{
    public:
    ViscosityKernel(Real r);
    Real W(Vec3d& rij)override;
    Vec3d GradW(Vec3d& r)override;
    Real LaplacianW(Vec3d& rij) override{
    Real r2 = SquaredNorm(rij);
    Real r1 = std::sqrt(r2);
    Real r3 = std::sqrt(r2)*r2;
    Real h6=h*h*h*h*h*h;
    return 45.0/(M_PI*h6)*(h-r1);
}

    protected:
    void SetVal()override;
};

