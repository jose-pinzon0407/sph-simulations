#pragma once
#define MORTON_BITS 7

using itype=uint32_t;

__constant__ int d_morton_bits;
struct cuMortonOrder{


    __device__
    static inline int MortonCode2d(int x, int y){
        // int z_val=0;

        // for(int i=0;i<MORTON_BITS;++i){
        //     z_val |= (int) ((x>>i) & 1) << (2*i);
        //     z_val |= (int) ((y>>i) & 1) << (2*i+1);
        // }

        // return z_val;

        return (Partition1By1(y)<<1)|Partition1By1(x);
    }

    __device__
    static inline itype MapThreeToOne(int c0, int c1, int c2, int nx, int ny, int nz){
        return MortonCode2d(c0,c1);
    }

    __device__
    static inline void MapOneToThree(int& c0, int& c1, int& c2, int id, int nx, int ny, int nz){
        // c0 = 0;
        // c1 = 0;
        // c2 = 0;

        // for(int i = 0; i < MORTON_BITS; ++i) {
        //     c0 |= ((id >> (2*i))     & 1) << i;
        //     c1 |= ((id >> (2*i + 1)) & 1) << i;
        // }
        c0 = Compact1By1(id >> 0);
        c1 = Compact1By1(id >> 1);
    }

    __device__ 
    static inline itype Partition1By1(itype x){
        x &= 0x0000ffff;

        x = (x ^ (x << 8)) & 0x00ff00ff;
        x = (x ^ (x << 4)) & 0x0f0f0f0f;
        x = (x ^ (x << 2)) & 0x33333333;
        x = (x ^ (x << 1)) & 0x55555555;
        return x;
    }

    __device__
    static inline itype Compact1By1(itype x){
        x &= 0x55555555;
        x = (x ^ (x>> 1)) & 0x33333333;
        x = (x ^ (x>> 2)) & 0x0f0f0f0f;
        x = (x ^ (x>> 4)) & 0x00ff00ff;
        x = (x ^ (x>> 8)) & 0x0000ffff;
        return x;

    }

};