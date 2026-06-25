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
    Real W(Vec3d& rij)override;
    Vec3d GradW(Vec3d& r)override;
    Real LaplacianW(Vec3d& r) override;

    protected:
    void SetVal()override;
};

class PressureKernel:public SmoothingFunction{
    public:
    PressureKernel(Real r);
    Real W(Vec3d& rij)override;
    Vec3d GradW(Vec3d& r)override;
    Real LaplacianW(Vec3d& r) override;

    protected:
    void SetVal()override;
};

class ViscosityKernel:public SmoothingFunction{
    public:
    ViscosityKernel(Real r);
    Real W(Vec3d& rij)override;
    Vec3d GradW(Vec3d& r)override;
    Real LaplacianW(Vec3d& r) override;

    protected:
    void SetVal()override;
};

