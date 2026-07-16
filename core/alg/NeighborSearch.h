#pragma once

#include "ConfigVariables.hpp"
#include "Domain.h"
#include "container/SortedContainer.h"
#include<fstream>
#include<functional>

struct NeighborSearch{

    public: 
    NeighborSearch(Domain& domain){
        SetGrid(domain);
    }
    void FindNeighbors(SortedSoAContainer& particles){
        size_t num_particles = particles.Size();

        for(size_t i=0;i<num_particles;++i){
            Real rxi = particles.GetParticles().rx[i];
            Real ryi = particles.GetParticles().ry[i];
            Real rzi = particles.GetParticles().rz[i];

            if(rxi < 0.0) std::cout<<"Out of domain\n";
            if(ryi < 0.0) std::cout<<"Out of domain\n";
            if(rzi < 0.0) std::cout<<"Out of domain\n";

            int ci0 = std::floor(rxi/grid.cell_size);
            int ci1 = std::floor(ryi/grid.cell_size);
            int ci2 = std::floor(rzi/grid.cell_size);

            ci0 = std::min(ci0, grid.nx-1);
            ci1 = std::min(ci1, grid.ny-1);

            // int ci0 = (rxi/grid.cell_size);
            // int ci1 = (ryi/grid.cell_size);
            // int ci2 = (rzi/grid.cell_size);

            uint idx = ci0 + ci1*grid.nx;
            particles.SetParticleCellIndex(i,idx);

        }

        // SortIntoCells

        particles.Build();
        OutputOrder(particles);
        //Create cell ranges
        SetCellRanges(num_particles,particles);
    }

    void OutputOrder(SortedSoAContainer& particles){
        std::ofstream output;            
        std::string name = "order_step_"+std::to_string(step)+".txt";
        output.open(name.c_str());
        for(size_t i=0;i<particles.Size();++i){
                output<<i<<"\t"<<particles.GetParticleIndex(i)<<"\t"<<particles.GetParticleCellIndex(i)<<"\t";
                output<<particles.pos(i)[0]<<"\t"<<particles.pos(i)[1];
                output<<std::endl;
        }
        output.close();
        name = "cell_order_step_"+std::to_string(step)+".txt";
        output.open(name.c_str());

        for(size_t i=0;i<grid.GetTotalCells();++i){
            output<<"Cell "<<i<<"["<<grid.cell_start[i]<<","<<grid.cell_end[i]<<"]"<<std::endl;
        }

        output.close();

        step++;
    }

    protected:
    int step=0;
    struct UniformGrid{
        Real cell_size;
        int nx,ny,nz;
        

        std::vector<int> cell_start, cell_end;

        unsigned int GetTotalCells(){
            return nx*ny*nz;
        }

        void ResizeCellRanges(){
            cell_start.resize(GetTotalCells(),-1);
            cell_end.resize(GetTotalCells(),-1);
        }

        void RefillCellRanges(){
            std::fill(cell_start.begin(),cell_start.end(),-1);
            std::fill(cell_end.begin(),cell_end.end(),-1);
        }

        int MapThreeToOne(const Vec3i& c){
            return c[0] +  c[1]*nx;// + c[2]*nx*ny;
        }

        Vec3i MapOneToThree(int cell_ind){
            int ci = cell_ind;
            Vec3i index;
            index[0]=ci%nx;
            index[1]=(ci/nx)%ny;
            index[2]=0;
            return index;
        }
    };

    UniformGrid grid;

    void SetGrid(Domain& domain){
        Vec3d box_length = domain.GetLength();

        grid.cell_size = H;
        grid.nx = std::ceil(box_length[0]/grid.cell_size);
        grid.ny = std::ceil(box_length[1]/grid.cell_size);
        grid.nz = std::max((int)std::ceil(box_length[2]/grid.cell_size),1);

        grid.ResizeCellRanges();

        std::cout<<"[NeighborSearch] Using uniform grid with: "<<grid.GetTotalCells()<<" cells."<<std::endl;
        std::cout<<"[NeighborSearch] Cells per direction: "<<grid.nx<<","<<grid.ny<<","<<grid.nz<<std::endl;
        std::cout<<"[NeighborSearch] Cell size: "<<grid.cell_size<<std::endl;
    }

    void SetCellRanges(uint N, SortedSoAContainer& particles){
        grid.RefillCellRanges();
        for(uint i=0;i<N;++i){
            uint cell_i_id = particles.GetParticleCellIndex(i);
            if(i==0 || cell_i_id != particles.GetParticleCellIndex(i-1)){
                grid.cell_start[cell_i_id] = i;
            }

            if(i==(N-1) || cell_i_id != particles.GetParticleCellIndex(i+1)){
                grid.cell_end[cell_i_id]= i+1;//open interval
            }
        }


        SortValidation(N);
    }

    void SortValidation(uint N){
        int accum_particles=0;
        for(size_t i=0;i<grid.GetTotalCells();++i){
            if(grid.cell_end[i] > -1 && grid.cell_start[i]> -1)
                accum_particles += (grid.cell_end[i]-grid.cell_start[i]);
        }
        std::cout<<"[NeighborSearch] Validation, counted particles "<<accum_particles<<" vs. N="<<N<<std::endl;
    }

    public: 
    int GetCellStart(int cell){
        return grid.cell_start[cell];
    }

    int GetCellEnd(int cell){
        return grid.cell_end[cell];
    }

    int Flatten(const Vec3i& c){
        return grid.MapThreeToOne(c);
    }

    Vec3i Unflatten(unsigned int index){
        return grid.MapOneToThree(index);
    }

    int GetCellsInDirectionX(){
        return grid.nx;
    }

    int GetCellsInDirectionY(){
        return grid.ny;
    }

    int GetCellsInDirectionZ(){
        return grid.nz;
    }

};