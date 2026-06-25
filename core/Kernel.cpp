#include "Kernel.h"


DensityKernel::DensityKernel(Real r):SmoothingFunction{r}{
    SetVal();
}

void DensityKernel::SetVal(){
    val=315.0/(64.0*M_PI*std::pow(h,9.0));
}

Real DensityKernel::W(Vec3d& rij){
    Real r2 = SquaredNorm(rij);
    return  (h2-r2)*(h2-r2)*(h2-r2)*val;
}

Vec3d DensityKernel::GradW(Vec3d& rij){
    Real r1 = VecNorm(rij);
    Real r2 = SquaredNorm(rij);
    Real dW = -6.0*val*(h*h-r2)*(h*h-r2);
    return (dW/r1)*rij;
}

Real DensityKernel::LaplacianW(Vec3d& rij){
    Real r2 = SquaredNorm(rij);
    return 6.0*val*(h*h-r2)*(7.0*r2-h*h);
}


/***************
 * PRESSURE KERNEL
 ***************/

PressureKernel::PressureKernel(Real r):SmoothingFunction{r}{
    SetVal();
}

void PressureKernel::SetVal(){
    val= 15.0/(M_PI*std::pow(h,6.0));
}

Real PressureKernel::W(Vec3d& rij){
    Real r1 = VecNorm(rij);
    return (h-r1)*(h-r1)*(h-r1)*val;
}

Vec3d PressureKernel::GradW(Vec3d& rij){
    Real r1 = VecNorm(rij);
    if(r1<TOL) return Vec3d{0.0,0.0,0.0};

    Real dW= -3.0*val*(h-r1)*(h-r1);
    return (dW/r1)*rij;
}

Real PressureKernel::LaplacianW(Vec3d& rij){
    //Not implemented
    return 0.0;
}

/***************
 * VISCOSITY KERNEL
 ***************/

ViscosityKernel::ViscosityKernel(Real r):SmoothingFunction{r}{
    SetVal();
}

void ViscosityKernel::SetVal(){
    val= 15.0/(2.0*M_PI*std::pow(h,3.0));
}

Real ViscosityKernel::W(Vec3d& rij){
    Real r1 = VecNorm(rij);
    Real w_value =-0.5*r1*r1*r1/(h*h*h)+r1*r1/(h*h)+0.5*h/r1-1;
    return  w_value*val;
}

Vec3d ViscosityKernel::GradW(Vec3d& rij){
    Real r2 = SquaredNorm(rij);
    Real r1 = std::sqrt(r2);
    Real dW= val*(-3.0*0.5*r2/(h*h*h)+2.0*r1/(h*h)-0.5*h/r2);
    return (dW/r1)*rij;
}

Real ViscosityKernel::LaplacianW(Vec3d& rij){
    //Not implemented
    return 0.0;
}