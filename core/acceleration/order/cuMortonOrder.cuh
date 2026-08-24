#pragma once
#define MORTON_BITS 7


__constant__ int d_morton_bits;
struct cuMortonOrder{


    __device__
    static inline int MortonCode2d(int x, int y){
        int z_val=0;

        for(int i=0;i<MORTON_BITS;++i){
            z_val |= (int) ((x>>i) & 1) << (2*i);
            z_val |= (int) ((y>>i) & 1) << (2*i+1);
        }

        return z_val;
    }

    __device__
    static inline int MapThreeToOne(int c0, int c1, int c2, int nx, int ny, int nz){
        return MortonCode2d(c0,c1);
    }

    __device__
    static void MapOneToThree(int& c0, int& c1, int& c2, int id, int nx, int ny, int nz){
        c0 = 0;
        c1 = 0;
        c2 = 0;

        for(int i = 0; i < MORTON_BITS; ++i) {
            c0 |= ((id >> (2*i))     & 1) << i;
            c1 |= ((id >> (2*i + 1)) & 1) << i;
        }
    }
};