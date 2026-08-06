#pragma once

struct MortonOrder{


    //Return xyz Morton codification
    int MortonCode(int x, int y, int z){
        int z_val=0;

        for(int i=0;i<21;++i){
        z_val |= (int)((x >> i) & 1) << (3 * i);
        z_val |= (int)((y >> i) & 1) << (3 * i + 1);
        z_val |= (int)((z >> i) & 1) << (3 * i + 2);
        }

        return z_val;
    }

    //Return xyz Morton codification
    int MortonCode2d(int x, int y){
        int z_val=0;

        for(int i=0;i<21;++i){
        z_val |= (int)((x >> i) & 1) << (2 * i);
        z_val |= (int)((y >> i) & 1) << (2 * i + 1);
        }

        return z_val;
    }

    //performs opposite to MortonCode to recover local indices
    void Deinterleave(int morton, int& x, int& y, int& z){
        
    }

};