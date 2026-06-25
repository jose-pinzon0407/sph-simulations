#include "Kernel.h"
#include "math/Math.h"
#include "Domain.h"
#include "container/Container.h"

#include "integrator/Integrator.h"
#include "container/Lattice.h"

int main(int argc, char** argv){
    
    Domain domain;
    domain.h_corner=Vec3d {10,20,0};
    domain.l_corner=Vec3d {0,0,0};
    int N=100;
    ParticleContainer container{N};
    ViscosityKernel kernel{H};
    
    Integrator integrator{0.01,&domain};

    LatticeCreator lattice;
    lattice.yVol=0.7;
    lattice.xVol=.5;
    lattice.dX=0.0625;
    lattice.dY=0.125;

    lattice.CreateLattice(domain,container);

    return 0;
}