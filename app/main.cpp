

// #include "compute/ListsForceComputation.h"
// #include "compute/GenericForceComputation.hpp"

// #include "alg/TypesNeighborSearch.h"
// #include "alg/GenericNeighborLists.hpp"

#include "integrator/Integrator.h"
#include "parallel/OMPFunctions.h"

// #include "utils/CheckPointWriter.h"
// #include "utils/TimerManager.h"
// #include "utils/InputParser.h"

#ifdef GPU
#include "acceleration/CUDAKernels.cuh"
#endif

// #include "algorithmic/MeasureAverageNeighborsPerParticle.h"
#include "simulation_configuration.h"
using namespace simulation;



// std::string StepHeader(int step, Container& particles, Domain& domain);


int main(int argc, char** argv){

    // InputParser parser;
    parser.ReadInputFile("example.txt");

    CreateTimers();

    #ifdef OPENMP
    OMPConfiguration();
    #endif

    #ifdef GPU
    std::cout<<"[Build] Using GPU acceleration\n";
    #endif

    ConfigureDomain();
    ConfigureLattice();
    ConfigureParameters();


    OutputOtherConfigurations();

    RunSimulation();



    return 0;
}
