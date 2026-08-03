#include "Kernel.h"
#include "Domain.h"

#include "container/MeasuringFunctions.h"
#include "container/Lattice.h"
#include "container/Container.h"
#include "container/SoAContainer.h"
#include "container/SortedContainer.h"

#include "alg/ForceComputation.h"
#include "alg/SoAForceComputation.h"
#include "alg/ListsForceComputation.h"
#include "alg/NeighborSearch.h"
#include "alg/NeighborLists.h"

#include "alg/EncodedListsForceComputation.h"

#include "integrator/Integrator.h"

#include "parallel/OMPFunctions.h"

#include "utils/CheckPointWriter.h"
#include "utils/TimerManager.h"
#include "utils/InputParser.h"

#ifdef GPU
#include "acceleration/CUDAKernels.cuh"
#endif

std::string StepHeader(int step, Container& particles, Domain& domain);


int main(int argc, char** argv){

    InputParser parser;
    parser.ReadInputFile("example.txt");

    TimerManager timers;

    timers.AddTimer("main_loop");
    timers.AddTimer("force_computation");
    timers.AddTimer("container_build");
    timers.Init();

    #ifdef OPENMP
    OMPConfiguration();
    #endif

    #ifdef GPU
    std::cout<<"[Build] Using GPU acceleration\n";
    #endif

    int steps=parser.GetIntParameter("steps");
    Real dt=parser.GetRealParameter("dt");
    Domain domain;
    {
        Real hx=parser.GetRealParameter("domain_x");
        Real hy=parser.GetRealParameter("domain_y");
        domain.h_corner=Vec3d{hx,hy,0};
        domain.l_corner=Vec3d{0,0,0};
    }


    NeighborSearch neighbor_search{domain,H};
    NeighborList list{H};
    // ParticleContainer container;
    SoAContainer container;

    ViscosityKernel kernel{H};
    
    Integrator integrator{dt,&domain};
    CheckPointWriter checkpoint{"positions_step_",parser.GetIntParameter("checkpoint_write_frequency")};
    LatticeCreator lattice;
    lattice.xVol=parser.GetRealParameter("lattice_xvol");
    lattice.yVol=parser.GetRealParameter("lattice_yvol");
    lattice.dX=parser.GetRealParameter("lattice_dx");
    lattice.dY=parser.GetRealParameter("lattice_dy");

    lattice.CreateLattice(domain,container);

    std::cout<<"[Configuration] Domain size: "<<domain.GetLength()<<"\n";
    std::cout<<"[Configuration] Particles="<<container.Size()<<"\n";
    std::cout<<"[Configuration] steps="<<steps<<"\n";
    std::cout<<"[Configuration] dt="<<dt<<"\n";
    std::cout<<"[Configuration] dX="<<lattice.dX<<"\n";
    std::cout<<"[Configuration] dY="<<lattice.dY<<"\n";
    std::cout<<"[Configuration] xVol="<<lattice.xVol<<"\n";
    std::cout<<"[Configuration] yVol="<<lattice.yVol<<"\n";
    std::cout<<"[Configuration] output frequency="<<parser.GetIntParameter("checkpoint_write_frequency")<<"\n";

    std::cout<<"[Configuration] Configuration complete."<<std::endl;

    neighbor_search.FindNeighbors(container);
    list.AssembleNeighborList(container,neighbor_search);


    std::cout<<"[Simulation] Output step 0 positions\n";
    checkpoint.WriteCheckPoint(0,container);
    

    // SoAForceComputation computation{H};
    // SortedForceComputation computation{H,neighbor_search};
    ListsForceComputation computation{H,neighbor_search,list};
    // EncodedListsForceComputation computation{H,neighbor_search,list};
    timers.Start("main_loop");
    for(int i=1;i<steps;++i){
        

        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);

        timers.Start("container_build");
        neighbor_search.FindNeighbors(container);
        list.AssembleNeighborList(container,neighbor_search);
        timers.Stop("container_build");

        checkpoint.WriteCheckPoint(i,container);


        timers.Start("force_computation");
        computation.ComputeForces(container);
        timers.Stop("force_computation");

        std::cout<<StepHeader(i,container,domain);

    }

    int N=container.Size();
    timers.Stop("main_loop");
    timers.OutputTimersToFile(N,steps);

    timers.PrintTimers();

    return 0;
}

std::string StepHeader(int step, Container& particles, Domain& domain){
    std::ostringstream out;

    out<<"Step= "<<step<<"\t N="<<ParticlesInDomain(domain,particles)<<"\t E_k="<<ComputeKineticEnergy(particles)<<"\t |v|_max="<<ComputeMaxVelocity(particles)<<std::endl;

    return out.str();
}