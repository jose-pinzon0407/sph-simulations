#include "Container.h"
#include "Domain.h"

struct LatticeCreator{

    double dX=0.1, dY=0.1, dZ=0.1;
    double xVol=1.0;
    double yVol=1.0;
    Vec3i particles_per_direction;
    int N;

    void CreateLattice(Domain& domain, Container& particles){
        particles_per_direction[0] = std::ceil(xVol*domain.GetLength()[0]/dX)+1;
        particles_per_direction[1] = (double)std::ceil(yVol*domain.GetLength()[1]/dY)+1;
        particles_per_direction[2] = std::max(domain.GetLength()[2]/dZ,1.0);

        N = (particles_per_direction[0])*(particles_per_direction[1])*particles_per_direction[2];

        particles.Resize(N);

        std::cout<<"[LatticeCreator] Total particles: "<<particles.Size()<<std::endl;
        std::cout<<"[LatticeCreator] Particles per direction: "<<particles_per_direction<<std::endl;

        int idx=0;
        for(int j=0;j<particles_per_direction[1];++j){
            for(int i=0;i<particles_per_direction[0];++i){
                idx = i + particles_per_direction[0]*j;
                particles.SetParticleIndex(idx,idx);
                particles.SetPosition(idx,dX*i,dY*j,0.0);
            }
        }

    }

};