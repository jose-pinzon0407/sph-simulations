#pragma once

#include<vector>
#include "container/SoAContainer.h"
#include "alg/NeighborSearch.h"
#include "encoders/Encoder.h"
#include "utils/Printer.h"

using neighbors = std::vector<std::vector<int>>;

/**
 * Essentially a lists tool
 */

struct AvgNeighborsPerParticle{
    
    Real median,deviation;
    int average, min=0, max=0, measured_steps=0;

    int num_cells=0, num_empty_cells=0;
    int avg_particles_per_cell=0.0, min_particles_per_cell=1e7, max_particles_per_cell=0;

    FileSystem file;

    const char* name = "neighborhood_and_cells.data";

    AvgNeighborsPerParticle(){
        bool created=file.CreateIfNotExists(name);
        if(created){
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"avg_npp"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"min_npp"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"max_npp"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"stddev_npp"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"#cells"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"#empty_cells"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"avg_ppc"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"min_ppc"<<"\t";
            file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<"max_ppc"<<"\n";
        }
        file.ReleaseStream();
    }


    void Measure(neighbors& list){
        int num_particles = list.size();
        int sum_neighbors=0;
        min=num_particles;

        for(int i=0;i<num_particles;++i){
            int i_neighbors= list[i].size();
            sum_neighbors += i_neighbors;
            max = std::max(max,i_neighbors);
            min = std::min(min,i_neighbors);
        }

        average = sum_neighbors/num_particles;

        for(int i=0;i<num_particles;++i){
            int i_neighbors= list[i].size();
            deviation += (i_neighbors - average)*(i_neighbors - average);
        }

        deviation = std::sqrt(deviation/num_particles);
        measured_steps++;

    }

    void MeasureCells(std::vector<int>& cell_start,std::vector<int>& cell_end){
        num_cells = cell_start.size();
        num_empty_cells=0;
        int sum_particles_per_cell=0;

        for(int i=0;i<num_cells;++i){
            int parts_in_cell = cell_end[i]-cell_start[i];
            sum_particles_per_cell += parts_in_cell;
            max_particles_per_cell = std::max(max_particles_per_cell,parts_in_cell);
            if(parts_in_cell>0){
                min_particles_per_cell = std::min(min_particles_per_cell,parts_in_cell);
            }
            
            if(parts_in_cell==0){
                num_empty_cells++;
            }
        }

        avg_particles_per_cell = sum_particles_per_cell/num_cells;

    }

    void PrintToFile(){

        file.OpenInputAppend(name);

        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<average<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<min<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<max<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<deviation<<"\t";

        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<num_cells<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<num_empty_cells<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<avg_particles_per_cell<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<min_particles_per_cell<<"\t";
        file.GetStream()<<std::fixed<<std::setw(20)<<std::left<<max_particles_per_cell<<"\n";

        file.ReleaseStream();
    }


};

