#pragma once

#include<vector>

/**
 * Delta codifier.
 * As a start we will implement differences between neighbors.
 */


struct Codifier{

    std::vector<std::vector<int>> neighbors_per_particle; 
    
    void Init(int N){
        neighbors_per_particle.resize(N);
    }

    void ResetCodifier(){
        int N = neighbors_per_particle.size();
        for(int i=0;i<N;++i){
            neighbors_per_particle[i].clear();
        }
    }

    //receives an ascending-sorted list of neigbor indices
    inline void Encode(int pi, std::vector<int>& particle_neighbors){
        int num_neighbors=particle_neighbors.size();
        
        std::vector<int>& i_deltas = neighbors_per_particle[pi];

        //store first neighbor verbatim
        i_deltas.emplace_back(particle_neighbors[0]);
        int first = particle_neighbors[0];

        //store differences
        for(int i=1;i<num_neighbors;++i){  
            int delta = particle_neighbors[i]-particle_neighbors[i-1]-1;
            // i_deltas.emplace_back(particle_neighbors[i] - first);
            i_deltas.emplace_back(delta);
        }


    }


};