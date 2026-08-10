#pragma once

#include<vector>
#include "container/SoAContainer.h"
#include "encoders/Encoder.h"

using neighbors = std::vector<std::vector<int>>;

template<class NeighborSearchType>
class GenericNeighborLists{

    public:
    GenericNeighborLists(double r);

    void AssembleNeighborList(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& neighbor_search);
    neighbors& GetNeighborList();
    std::vector<int>& GetNeighborsOfParticle(int i){
        return neighbor_lists[i];
    }
    void ResetNeighborList();
    void Resize(int N);
    
    void SeeList(){
        for(int i=0;i<neighbor_lists.size();++i){
            std::cout<<"Particle "<<i<<" has neighbors "<<neighbor_lists[i].size()<<": ";
            for(int j=0;j<neighbor_lists[i].size();++j){
                std::cout<<neighbor_lists[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
    }

    Codifier& GetCodifier(){
        return encoder;
    }

    protected:
    Real h2;
    std::vector<std::vector<int>> neighbor_lists;
    Codifier encoder;

};

#include "GenericNeighborLists.tpp"