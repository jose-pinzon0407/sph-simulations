#pragma once

#include "ConfigVariables.hpp"
#include "Domain.h"
#include "container/SortedContainer.h"
#include<fstream>
#include<functional>

class NeighborSearch{

    public: 
    NeighborSearch(Domain& domain, Real cell_size);
    void FindNeighbors(SoAContainer& particles);

    std::vector<int>& GetEnd();
    std::vector<int>& GetStart();
    int GetCellStart(int cell);
    int GetCellEnd(int cell);
    int Flatten(const Vec3i& c);
    Vec3i Unflatten(unsigned int index);
    int GetCellsInDirectionX();
    int GetCellsInDirectionY();
    int GetCellsInDirectionZ();
    int GetTotalCells();

    protected:

    void OutputOrder(SoAContainer& particles);
    void SetGrid(Domain& domain, Real radius);
    void SetCellRanges(uint N, SoAContainer& particles);
    void SortValidation(uint N);
    void SortParticleIds(SoAContainer& container);
    void SortContainer(SoAContainer& container);
    void Sort(SoAContainer& container);
    void SetParticleCellIds(SoAContainer& container);

    private:


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
    

};