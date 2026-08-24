#pragma once

/**
 * Functions that configure the simulation in the main loop
 */

#include "integrator/Integrator.h"

#include "algorithmic/MeasureAverageNeighborsPerParticle.h"

#include "utils/TimerManager.h"
#include "utils/InputParser.h"
#include "utils/CheckPointWriter.h"

#include "Kernel.h"
#include "Domain.h"

#include "container/MeasuringFunctions.h"
#include "container/Lattice.h"
#include "container/Container.h"
#include "container/SoAContainer.h"

#include "compute/ListsForceComputation.h"

#include "alg/TypesNeighborSearch.h"
#include "compute/GenericForceComputation.hpp"
#include "alg/GenericNeighborLists.hpp"

#include "algorithmic/MultipassNeighborLists.h"


namespace simulation{

    void PrintCellsToFile(int step, std::vector<int>& start, std::vector<int>& end){
    std::ofstream file;
    std::string name="cells_"+std::to_string(step)+".txt";
    file.open(name.c_str());

    for(int i=0;i<start.size();++i){
        file<<"cell:"<<i<<"["<<start[i]<<","<<end[i]<<")\n";
    }

    file.close();

    }


    void PrintCellIDToFile(int step, std::vector<int>& cell_id){
    std::ofstream file;
    std::string name="ids_"+std::to_string(step)+".txt";
    file.open(name.c_str());

    for(int i=0;i<cell_id.size();++i){
        file<<"particle:"<<i<<" in cell: "<<cell_id[i]<<"\n";
    }

    file.close();

    }

InputParser parser;
Domain domain;
TimerManager timers;
SoAContainer container;
int N, steps;
Real dt;
AvgNeighborsPerParticle measure;
CheckPointWriter* checkpoint=nullptr;

void CreateTimers(){

    timers.AddTimer("main_loop");
    timers.AddTimer("force_computation");
    timers.AddTimer("calculate_pressure");
    timers.AddTimer("approximate_density");
    timers.AddTimer("neighbor_search");
    timers.AddTimer("neighbor_list");
    timers.AddTimer("integration");
    timers.AddTimer("IO");

    #ifdef GPU
    timers.AddTimer("cuda_buffer_copy");
    timers.AddTimer("buffer_inversion");
    #endif


    timers.Init();

}

void ConfigureDomain(){
    Real hx=parser.GetRealParameter("domain_x");
    Real hy=parser.GetRealParameter("domain_y");
    domain.h_corner=Vec3d{hx,hy,0};
    domain.l_corner=Vec3d{0,0,0};

    
}

void ConfigureLattice(){
    LatticeCreator lattice;
    lattice.xVol=parser.GetRealParameter("lattice_xvol");
    lattice.yVol=parser.GetRealParameter("lattice_yvol");
    lattice.dX=parser.GetRealParameter("lattice_dx");
    lattice.dY=parser.GetRealParameter("lattice_dy");

    std::cout<<"[Configuration] dX="<<lattice.dX<<"\n";
    std::cout<<"[Configuration] dY="<<lattice.dY<<"\n";
    std::cout<<"[Configuration] xVol="<<lattice.xVol<<"\n";
    std::cout<<"[Configuration] yVol="<<lattice.yVol<<"\n";
    
    lattice.CreateLattice(domain,container);
    std::cout<<"[Configuration] Lattice created in domain\n";
}

void ConfigureParameters(){
    steps=parser.GetIntParameter("steps");
    dt=parser.GetRealParameter("dt");

    checkpoint = new CheckPointWriter{"positions_step_",parser.GetIntParameter("checkpoint_write_frequency")};
}

void OutputPerformanceProfile(int steps, int N){
    timers.PrintTimers();
    timers.OutputTimersToFile(N,steps);
}


void OutputOtherConfigurations(){
    std::cout<<"[Configuration] Domain size: "<<domain.GetLength()<<"\n";
    std::cout<<"[Configuration] Particles="<<container.Size()<<"\n";
    std::cout<<"[Configuration] steps="<<steps<<"\n";
    std::cout<<"[Configuration] dt="<<dt<<"\n";
    std::cout<<"[Configuration] output frequency="<<parser.GetIntParameter("checkpoint_write_frequency")<<"\n";

    std::cout<<"[Configuration] Configuration complete."<<std::endl;
}

std::string StepHeader(int step, Container& particles, Domain& domain){
    std::ostringstream out;

    out<<"Step= "<<step<<"\t N="<<ParticlesInDomain(domain,particles)<<"\t E_k="<<ComputeKineticEnergy(particles)<<"\t |v|_max="<<ComputeMaxVelocity(particles)<<std::endl;

    return out.str();
}

void End(){
    int N=container.Size();
    OutputPerformanceProfile(steps, N);
    measure.PrintToFile();
}


template<class NeighborSearchType>
void PreLoopComputation(GenericNeighborSearch<NeighborSearchType>& neighbor_search, GenericNeighborLists<NeighborSearchType>& neighbor_list, GenericForceComputation<NeighborSearchType>& computation ){

    timers.Start("neighbor_search");
    neighbor_search.FindNeighbors(container);
    timers.Stop("neighbor_search");

    // PrintCellsToFile(0,neighbor_search.GetStart(),neighbor_search.GetEnd());

    // timers.Start("neighbor_list");
    // neighbor_list.AssembleNeighborList(container,neighbor_search);
    // timers.Stop("neighbor_list");


    std::cout<<"[Simulation] Output step 0 positions\n";
    // checkpoint->WriteCheckPoint(0,container);

    timers.Start("approximate_density");
    computation.ApproximateDensity(container);
    timers.Stop("approximate_density");

    timers.Start("calculate_pressure");
    computation.CalculatePressure(container);
    timers.Stop("calculate_pressure");

    timers.Start("force_computation");
    computation.ComputeForces(container);
    timers.Stop("force_computation");

    timers.Start("IO");
    std::cout<<StepHeader(0,container,domain);
    timers.Stop("IO");

}

template<class NeighborSearchType>
void PreLoopComputationWithLists(GenericNeighborSearch<NeighborSearchType> & neighbor_search, GenericNeighborLists<NeighborSearchType> & neighbor_list, ListsForceComputation& computation ){
    timers.Start("main_loop");
    timers.Start("neighbor_search");
    neighbor_search.FindNeighbors(container);
    timers.Stop("neighbor_search");

    timers.Start("neighbor_list");
    neighbor_list.AssembleNeighborList(container,neighbor_search);
    timers.Stop("neighbor_list");


    std::cout<<"[Simulation] Output step 0 positions\n";
    // checkpoint->WriteCheckPoint(0,container);

    timers.Start("approximate_density");
    computation.ApproximateDensity(container);
    timers.Stop("approximate_density");

    timers.Start("calculate_pressure");
    computation.CalculatePressure(container);
    timers.Stop("calculate_pressure");

    timers.Start("force_computation");
    computation.ComputeForces(container);
    timers.Stop("force_computation");

    timers.Start("IO");
    std::cout<<StepHeader(0,container,domain);
    timers.Stop("IO");
    timers.Stop("main_loop");


}

template<class NeighborSearchType>
void MainLoopComputationWithLists(GenericNeighborSearch<NeighborSearchType>& neighbor_search, GenericNeighborLists<NeighborSearchType>& neighbor_list, ListsForceComputation& computation){
    Integrator integrator{dt,&domain};

    timers.Start("main_loop");
    for(int i=1;i<steps;++i){
        
        timers.Start("integration");
        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);
        timers.Stop("integration");


        timers.Start("neighbor_search");
        neighbor_search.FindNeighbors(container);
        timers.Stop("neighbor_search");

        timers.Start("neighbor_list");
        neighbor_list.AssembleNeighborList(container,neighbor_search);
        timers.Stop("neighbor_list");

        timers.Start("IO");
        checkpoint->WriteCheckPoint(i,container);
        timers.Stop("IO");

        timers.Start("approximate_density");
        computation.ApproximateDensity(container);
        timers.Stop("approximate_density");

        timers.Start("calculate_pressure");
        computation.CalculatePressure(container);
        timers.Stop("calculate_pressure");

        timers.Start("force_computation");
        computation.ComputeForces(container);
        timers.Stop("force_computation");

        // measure.Measure(neighbor_list.GetNeighborList());
        // measure.MeasureCells(neighbor_search.GetStart(),neighbor_search.GetEnd());

        timers.Start("IO");
        std::cout<<StepHeader(i,container,domain);
        timers.Stop("IO");

    }

    timers.Stop("main_loop");

}

template<class NeighborSearchType>
void PreLoopComputationWithMPLists(GenericNeighborSearch<NeighborSearchType> & neighbor_search, MultipassNeighborLists<NeighborSearchType> & neighbor_list, MPListsForceComputation& computation ){

    timers.Start("neighbor_search");
    neighbor_search.FindNeighbors(container);
    timers.Stop("neighbor_search");

    timers.Start("neighbor_list");
    neighbor_list.AssembleNeighborList(container,neighbor_search);
    timers.Stop("neighbor_list");


    std::cout<<"[Simulation] Output step 0 positions\n";
    // checkpoint->WriteCheckPoint(0,container);

    timers.Start("approximate_density");
    computation.ApproximateDensity(container);
    timers.Stop("approximate_density");

    timers.Start("calculate_pressure");
    computation.CalculatePressure(container);
    timers.Stop("calculate_pressure");

    timers.Start("force_computation");
    computation.ComputeForces(container);
    timers.Stop("force_computation");

    timers.Start("IO");
    std::cout<<StepHeader(0,container,domain);
    timers.Stop("IO");

}

template<class NeighborSearchType>
void MainLoopComputationWithMPLists(GenericNeighborSearch<NeighborSearchType>& neighbor_search, MultipassNeighborLists<NeighborSearchType>& neighbor_list, MPListsForceComputation& computation){
    Integrator integrator{dt,&domain};

    timers.Start("main_loop");
    for(int i=1;i<steps;++i){
        
        timers.Start("integration");
        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);
        timers.Stop("integration");


        timers.Start("neighbor_search");
        neighbor_search.FindNeighbors(container);
        timers.Stop("neighbor_search");

        timers.Start("neighbor_list");
        neighbor_list.AssembleNeighborList(container,neighbor_search);
        timers.Stop("neighbor_list");

        timers.Start("IO");
        checkpoint->WriteCheckPoint(i,container);
        timers.Stop("IO");

        timers.Start("approximate_density");
        computation.ApproximateDensity(container);
        timers.Stop("approximate_density");

        timers.Start("calculate_pressure");
        computation.CalculatePressure(container);
        timers.Stop("calculate_pressure");

        timers.Start("force_computation");
        computation.ComputeForces(container);
        timers.Stop("force_computation");

        // measure.Measure(neighbor_list.GetNeighborList(),neighbor_list.GetOffsetPerParticle());
        // measure.MeasureCells(neighbor_search.GetStart(),neighbor_search.GetEnd());

        timers.Start("IO");
        std::cout<<StepHeader(i,container,domain);
        timers.Stop("IO");

    }

    timers.Stop("main_loop");

}

template<class NeighborSearchType>
void MainLoopComputation(GenericNeighborSearch<NeighborSearchType>& neighbor_search, GenericNeighborLists<NeighborSearchType>& neighbor_list, GenericForceComputation<NeighborSearchType>& computation){
    Integrator integrator{dt,&domain};

    timers.Start("main_loop");
    for(int i=1;i<steps;++i){
        
        timers.Start("integration");
        integrator.UpdateVelocities(container);
        integrator.UpdatePositions(container);
        timers.Stop("integration");


        timers.Start("neighbor_search");
        neighbor_search.FindNeighbors(container);
        timers.Stop("neighbor_search");

        // PrintCellsToFile(i,neighbor_search.GetStart(),neighbor_search.GetEnd());


        // timers.Start("neighbor_list");
        // neighbor_list.AssembleNeighborList(container,neighbor_search);
        // timers.Stop("neighbor_list");

        timers.Start("IO");
        checkpoint->WriteCheckPoint(i,container);
        timers.Stop("IO");

        timers.Start("approximate_density");
        computation.ApproximateDensity(container);
        timers.Stop("approximate_density");

        timers.Start("calculate_pressure");
        computation.CalculatePressure(container);
        timers.Stop("calculate_pressure");

        timers.Start("force_computation");
        computation.ComputeForces(container);
        timers.Stop("force_computation");

        // measure.Measure(neighbor_list.GetNeighborList());
        // measure.MeasureCells(neighbor_search.GetStart(),neighbor_search.GetEnd());

        timers.Start("IO");
        std::cout<<StepHeader(i,container,domain);
        timers.Stop("IO");

    }

    timers.Stop("main_loop");

}

template<class NeighborSearchType>
void RunSimulationWithCellNeighbors(){

    GenericNeighborSearch<NeighborSearchType> gns{domain,H};
    GenericNeighborLists<NeighborSearchType> list{H};
    GenericForceComputation<NeighborSearchType> computation{H,gns};

    
    PreLoopComputation<NeighborSearchType>(gns,list,computation);

    MainLoopComputation<NeighborSearchType>(gns,list,computation);

    // End();

}

template<class NeighborSearchType>
void RunSimulationWithLists(){

    GenericNeighborSearch<NeighborSearchType> gns{domain,H};
    GenericNeighborLists<NeighborSearchType> list{H};
    ListsForceComputation computation{H,list.GetNeighborList()};

    
    PreLoopComputationWithLists<NeighborSearchType>(gns,list,computation);

    MainLoopComputationWithLists<NeighborSearchType>(gns,list,computation);

    // End();

}

template<class NeighborSearchType>
void RunSimulationWithMPLists(){
    GenericNeighborSearch<NeighborSearchType> gns{domain,H};
    MultipassNeighborLists<NeighborSearchType> list{H};
    MPListsForceComputation computation{H,list.GetNeighborList(),list.GetOffsetPerParticle()};

    
    PreLoopComputationWithMPLists<NeighborSearchType>(gns,list,computation);

    MainLoopComputationWithMPLists<NeighborSearchType>(gns,list,computation);

    // End();
}



//Launcher function
void RunSimulation(){


    if(parser.GetParameter("compute_type")=="mp-lists"){
        if(parser.GetParameter("order")=="morton"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithMPLists<MortonOrder>();
            return;
        }


        if(parser.GetParameter("order")=="rowmajor"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithMPLists<RowMajorOrder>();
            return;
            
        }
    }

    if(parser.GetParameter("compute_type")=="sorted"){
        std::cout<<"[Configuration] Using "<<parser.GetParameter("compute_type")<<" computation type\n";

        if(parser.GetParameter("order")=="morton"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithCellNeighbors<MortonOrder>();
            return;
        }


        if(parser.GetParameter("order")=="rowmajor"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithCellNeighbors<RowMajorOrder>();
            
            return;
            
        }

    }
    else{
        std::cout<<"[Configuration] Using "<<parser.GetParameter("compute_type")<<" computation type\n";

        if(parser.GetParameter("order")=="morton"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithLists<MortonOrder>();
            return;
        }


        if(parser.GetParameter("order")=="rowmajor"){
            std::cout<<"[Configuration] Using "<<parser.GetParameter("order")<<" ordering\n";
            RunSimulationWithLists<RowMajorOrder>();
            return;
            
        }
    }


}


}//end simulation namespace


