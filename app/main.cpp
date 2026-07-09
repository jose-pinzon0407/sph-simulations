#include "Kernel.h"
#include "math/Math.h"
#include "Domain.h"
#include "container/Container.h"

#include "integrator/Integrator.h"
#include "container/Lattice.h"

#include "utils/CheckPointWriter.h"
#include "utils/Timer.h"
#include "alg/ForceComputation.h"
#include "container/MeasuringFunctions.h"

std::string StepHeader(int step, Container& particles, Domain& domain);

Timer force_computation{"force_computation"};
Timer density_approximation{"density_approximation"};
Timer main_loop{"main_loop"};

int main(int argc, char** argv){
    int steps=100;
    Real dt=0.005;
    Domain domain;
    domain.h_corner=Vec3d {10,20,0};
    domain.l_corner=Vec3d {0,0,0};
    int N=10;
    ParticleContainer container;
    ViscosityKernel kernel{H};
    
    Integrator integrator{dt,&domain};

    LatticeCreator lattice;
    lattice.yVol=0.7;
    lattice.xVol=.5;
    lattice.dX=0.25;
    lattice.dY=0.25;

    lattice.CreateLattice(domain,container);

    CheckPointWriter checkpoint{"positions_step_",1};
    checkpoint.WriteCheckPoint(0,container);

    ForceComputation computation{H,container.Size()};
    // computation.ApproximateDensity(container);
    // computation.CalculatePressure(container);
    computation.ComputeForces(container);

    main_loop.Start();
    for(int i=0;i<steps;++i){
        
        std::cout<<StepHeader(i,container,domain);

        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);
        checkpoint.WriteCheckPoint(i,container);

        // density_approximation.Start();
        // computation.ApproximateDensity(container);
        // density_approximation.Stop();
// 
        // computation.CalculatePressure(container);

        force_computation.Start();
        computation.ComputeForces(container);
        force_computation.Stop();

    }
    main_loop.Stop();

    std::cout<<"Main loop ="<<main_loop.GetTotalTimeInSeconds()<<" s"<<std::endl;
    std::cout<<force_computation.GetName()<<"="<<force_computation.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<density_approximation.GetName()<<"="<<density_approximation.GetTotalTimeInSeconds()<<"s"<<std::endl;

    return 0;
}

std::string StepHeader(int step, Container& particles, Domain& domain){
    std::ostringstream out;

    out<<"Step= "<<step<<"\t N="<<ParticlesInDomain(domain,particles)<<"\t E_k="<<ComputeKineticEnergy(particles)<<"\t |v|_max="<<ComputeMaxVelocity(particles)<<std::endl;

    return out.str();
}