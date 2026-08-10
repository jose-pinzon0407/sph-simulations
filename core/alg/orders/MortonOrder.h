#pragma once

#include "ConfigVariables.hpp"
#include "Domain.h"

struct MortonOrder{

    Real cell_size;
    int nx, ny, nz;

    std::vector<int> cell_start, cell_end;

    std::vector<int> cell_x, cell_y, cell_z;

    void Initialize(Domain& domain, Real radius){
        cell_size=radius;
        Vec3d box_length=domain.GetLength();

        nx = std::ceil(box_length[0]/cell_size);
        ny = std::ceil(box_length[1]/cell_size);
        nz = std::max((int)std::ceil(box_length[2]/cell_size),1);

        AddCellPadding();
        int num_cells = GetTotalCells();

        cell_x.resize(num_cells);
        cell_y.resize(num_cells);
        cell_z.resize(num_cells);

        for(int cz=0;cz<nz;++cz){
            for(int cy=0;cy<ny;++cy){
                for(int cx=0;cx<nx;++cx){
                    int morton_code = MortonCode2d(cx,cy);
                    cell_x[morton_code]=cx;
                    cell_y[morton_code]=cy;
                    cell_z[morton_code]=cz;
            }
        }
    }

        ResizeCellRanges();

        std::cout<<"[MortonOrder] Using uniform grid with: "<<GetTotalCells()<<" cells."<<std::endl;
        std::cout<<"[MortonOrder] Cells per direction: "<<nx<<","<<ny<<","<<nz<<std::endl;
        std::cout<<"[MortonOrder] Cell size: "<<cell_size<<std::endl;

    }

    unsigned int GetTotalCells(){
        return nx*ny*nz;
    }

    void ResizeCellRanges(){
        std::cout<<"Resize to= "<<GetTotalCells()<<std::endl;
        int num_cells=GetTotalCells();
        cell_start.resize(num_cells,-1);
        cell_end.resize(num_cells,-1);
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
        return MortonCode2d(x,y);
    }

    Vec3i MapOneToThree(int cell_ind){
        Vec3i index;
        index[0]=cell_x[cell_ind];
        index[1]=cell_y[cell_ind];
        index[2]=cell_z[cell_ind];
        return index;
    }

    int NexPowerOfTwo(int n){
        int power=1;

        while(power <n ){
            power *= 2;
        }

        return power;
    }

    void AddCellPadding(){
        std::cout<<"[MortonOrder] Padding of the total cells is required\n";
        int old_num_cells = GetTotalCells();
        std::cout<<"[MortonOrder] Given configuration is: ("<<nx<<","<<ny<<","<<nz<<") a total of: "<<old_num_cells<<" cells\n"; 
        
        nx = NexPowerOfTwo(nx);
        ny = NexPowerOfTwo(ny);
        nz = NexPowerOfTwo(nz);
        std::cout<<"[MortonOrder] After padding the configuration is: ("<<nx<<","<<ny<<","<<nz<<") a total of: "<<GetTotalCells()<<" cells\n";
        std::cout<<"[MortonOrder] Added empty cells for padding: "<<GetTotalCells()-old_num_cells<<"\n";

    }

    int MortonCode3d(int x, int y, int z){
        int z_val=0;

        for(int i=0;i<21;++i){
        z_val |= (int)((x >> i) & 1) << (3 * i);
        z_val |= (int)((y >> i) & 1) << (3 * i + 1);
        z_val |= (int)((z >> i) & 1) << (3 * i + 2);
        }

        return z_val;
    }

    int MortonCode2d(int x, int y){
        int z_val=0;

        for(int i=0;i<21;++i){
        z_val |= (int)((x >> i) & 1) << (2 * i);
        z_val |= (int)((y >> i) & 1) << (2 * i + 1);
        }

        return z_val;
    }

};
