#pragma once

#include<vector>
#include "container/SoAContainer.h"


template<class NeighborSearchType>
class MultipassNeighborLists{

    public: 
    MultipassNeighborLists(Real r):h2{r*r}{

    }
    void AssembleNeighborList(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& ns){
        CountNeighborsPerParticle(container,ns);
        CalculateOffsets();
        AssignIndices(container,ns);
    }

    void ResetNeighborList(){


    }
    void Resize(int N){
        neighbor_list.resize(N);
    }

    std::vector<int>& GetOffsetPerParticle(){
        return offset_per_particle;
    }

    int GetOffsetFromParticle(int i){
        return offset_per_particle[i];
    }

    std::vector<int>& GetNeighborList(){
        return neighbor_list;
    }


    protected:

    void StoreNeighborIndicesPerParticle(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& ns){

    }

    void CalculateOffsets(){

        size_t N=num_neighbors_per_particle.size();
        offset_per_particle.resize(N+1);
        total_neighbors=0;
        offset_per_particle[0] = 0;
        
        for(int i=0;i<N;++i){
            offset_per_particle[i+1] = num_neighbors_per_particle[i] + offset_per_particle[i];
        }

        total_neighbors =  offset_per_particle[N];
        Resize(total_neighbors);

    }

    void CountNeighborsPerParticle(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& ns){
        size_t N = container.Size();
        num_neighbors_per_particle.resize(N);

        #pragma omp parallel for schedule(runtime)
        for(size_t i=0;i<N;++i){
            int count_i=0;
            Real rxi = container.GetParticles().rx[i];
            Real ryi = container.GetParticles().ry[i];
            Real rzi = container.GetParticles().rz[i];

            uint cell_i = container.GetParticleCellIndex(i);
            Vec3i m_idx_i = ns.Unflatten(cell_i);

            for(int c1=-1;c1<=1;++c1){
                for(int c2=-1;c2<=1;++c2){
                
                Vec3i m_idx_j = m_idx_i;
                
                m_idx_j[0] += c1;
                
                if(m_idx_j[0] < 0 || m_idx_j[0] >= ns.GetCellsInDirectionX()) continue;
                
                m_idx_j[1] += c2;
                
                if(m_idx_j[1] < 0|| m_idx_j[1] >= ns.GetCellsInDirectionY()) continue;
                

                int cell_j = ns.Flatten(m_idx_j);


                int j_start = ns.GetCellStart(cell_j);
                int j_end = ns.GetCellEnd(cell_j);

                if(j_start < 0 ||  j_end<0){
                    continue;
                }
                
                for(int pj=j_start;pj<j_end;++pj){

                    Real dx = rxi-container.GetParticles().rx[pj];
                    Real dy = ryi-container.GetParticles().ry[pj];
                    Real dz = rzi-container.GetParticles().rz[pj];

                    Real r2 = dx*dx+dy*dy+dz*dz; 

                    if(r2>h2) continue;

                    count_i++;


                }
                }   
            }

            num_neighbors_per_particle[i]=count_i;
        }

    }

    void AssignIndices(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& ns){

        int N=container.Size();

        #pragma omp parallel for schedule(runtime)
        for(size_t i=0;i<N;++i){
            Real rxi = container.GetParticles().rx[i];
            Real ryi = container.GetParticles().ry[i];
            Real rzi = container.GetParticles().rz[i];

            int offset_i = offset_per_particle[i];

            uint cell_i = container.GetParticleCellIndex(i);
            Vec3i m_idx_i = ns.Unflatten(cell_i);

            for(int c1=-1;c1<=1;++c1){
                for(int c2=-1;c2<=1;++c2){
                
                    Vec3i m_idx_j = m_idx_i;
                    
                    m_idx_j[0] += c1;
                    
                    if(m_idx_j[0] < 0 || m_idx_j[0] >= ns.GetCellsInDirectionX()) continue;
                    
                    m_idx_j[1] += c2;
                    
                    if(m_idx_j[1] < 0|| m_idx_j[1] >= ns.GetCellsInDirectionY()) continue;
                    

                    int cell_j = ns.Flatten(m_idx_j);


                    int j_start = ns.GetCellStart(cell_j);
                    int j_end = ns.GetCellEnd(cell_j);

                    if(j_start < 0 ||  j_end<0){
                        continue;
                    }
                    
                    for(int pj=j_start;pj<j_end;++pj){

                        Real dx = rxi-container.GetParticles().rx[pj];
                        Real dy = ryi-container.GetParticles().ry[pj];
                        Real dz = rzi-container.GetParticles().rz[pj];

                        Real r2 = dx*dx+dy*dy+dz*dz; 

                        if(r2>h2) continue;

                        neighbor_list[offset_i] = pj;
                        offset_i++;


                    }
                }   
            }

        }
    }

    Real h2;
    int total_neighbors=0;
    std::vector<int> num_neighbors_per_particle;
    std::vector<int> offset_per_particle;
    std::vector<int> neighbor_list;
    
};