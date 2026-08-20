#include "integrator/Integrator.h"
#include "parallel/OMPFunctions.h"


#ifdef GPU
#include "cu_simulation_configuration.h"
#else
#include "simulation_configuration.h"
#endif

using namespace simulation;
#ifdef GPU
using namespace gpu_simulation;
#endif


// std::string StepHeader(int step, Container& particles, Domain& domain);


int main(int argc, char** argv){
    std::cout<<"[Smoothed Particle Hydrodynamics]\n";
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

    #ifdef GPU
    cudaConfigure(domain,H,parser.GetIntParameter("gpu_threads_per_block"),container.Size());
    cudaSimulationSetup(container);
    #endif

    OutputOtherConfigurations();

    #ifdef GPU
    RunSimulationGPU();
    #else
    RunSimulation();
    #endif

    End();


    return 0;
}
