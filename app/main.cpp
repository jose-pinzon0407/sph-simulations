#include "Kernel.h"
#include "Domain.h"

#include "container/MeasuringFunctions.h"
#include "container/Lattice.h"
#include "container/Container.h"
#include "container/SoAContainer.h"
#include "container/SortedContainer.h"

#include "alg/ForceComputation.h"
#include "alg/SoAForceComputation.h"
#include "alg/NeighborSearch.h"

#include "integrator/Integrator.h"


#include "utils/CheckPointWriter.h"
#include "utils/Timer.h"

std::string StepHeader(int step, Container& particles, Domain& domain);

Timer force_computation{"force_computation"};
Timer density_approximation{"density_approximation"};
Timer container_build{"container_build"};
Timer main_loop{"main_loop"};

int main(int argc, char** argv){
    int steps=5;
    Real dt=0.005;
    Domain domain;
    domain.h_corner=Vec3d {10,20,0};
    domain.l_corner=Vec3d {0,0,0};
    NeighborSearch neighbor_search{domain,H};

    std::cout<<"[Domain] Length:"<<domain.GetLength()<<std::endl;

    int N=10;
    // ParticleContainer container;
    // SoAContainer container;
    SortedSoAContainer container;
    ViscosityKernel kernel{H};
    
    Integrator integrator{dt,&domain};

    LatticeCreator lattice;
    lattice.yVol=0.9;
    lattice.xVol=.5;
    lattice.dX=0.25;
    lattice.dY=0.25;

    lattice.CreateLattice(domain,container);
    neighbor_search.FindNeighbors(container);

    CheckPointWriter checkpoint{"positions_step_",0};
    checkpoint.WriteCheckPoint(0,container);


    // SoAForceComputation computation{H};
    SortedForceComputation computation{H,neighbor_search};


    main_loop.Start();
    for(int i=1;i<steps;++i){
        
        std::cout<<StepHeader(i,container,domain);

        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);

        container_build.Start();
        neighbor_search.FindNeighbors(container);
        container_build.Stop();

        checkpoint.WriteCheckPoint(i,container);
    

        force_computation.Start();
        computation.ComputeForces(container);
        force_computation.Stop();

    }
    main_loop.Stop();

    std::cout<<"Main loop ="<<main_loop.GetTotalTimeInSeconds()<<" s"<<std::endl;
    std::cout<<force_computation.GetName()<<"="<<force_computation.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<density_approximation.GetName()<<"="<<density_approximation.GetTotalTimeInSeconds()<<"s"<<std::endl;
    std::cout<<container_build.GetName()<<"="<<container_build.GetTotalTimeInSeconds()<<"s"<<std::endl;

    return 0;
}

std::string StepHeader(int step, Container& particles, Domain& domain){
    std::ostringstream out;

    out<<"Step= "<<step<<"\t N="<<ParticlesInDomain(domain,particles)<<"\t E_k="<<ComputeKineticEnergy(particles)<<"\t |v|_max="<<ComputeMaxVelocity(particles)<<std::endl;

    return out.str();
}