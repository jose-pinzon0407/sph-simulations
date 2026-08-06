#pragma once

#include "ConfigVariables.hpp"
#include "Domain.h"

struct RowMajorOrder{

    Real cell_size;
    int nx, ny, nz;

    std::vector<int> cell_start, cell_end;

    void Initialize(Domain& domain, Real radius){
        cell_size=radius;
        Vec3d box_length=domain.GetLength();

        nx = std::ceil(box_length[0]/cell_size);
        ny = std::ceil(box_length[1]/cell_size);
        nz = std::max((int)std::ceil(box_length[2]/cell_size),1);

        ResizeCellRanges();

        std::cout<<"[RowMajorOrder] Using uniform grid with: "<<GetTotalCells()<<" cells."<<std::endl;
        std::cout<<"[RowMajorOrder] Cells per direction: "<<nx<<","<<ny<<","<<nz<<std::endl;
        std::cout<<"[RowMajorOrder] Cell size: "<<cell_size<<std::endl;

    }

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
        return MapThreeToOne(c[0],c[1],c[2]);
        // return c[0] +  c[1]*nx;// + c[2]*nx*ny;
    }

    int MapThreeToOne(int x, int y, int z){
        return x + y*nx;
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
