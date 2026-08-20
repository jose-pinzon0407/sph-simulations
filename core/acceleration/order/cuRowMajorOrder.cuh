#pragma once

struct cuRowMajorOrder{
    __device__
    static inline int MapThreeToOne(int c0, int c1, int c2, int nx, int ny, int nz){
        return c0 + c1*nx;
    }

    __device__
    static void MapOneToThree(int& c0, int& c1, int& c2, int id, int nx, int ny, int nz){
        c0 = id%nx;
        c1 = (id/nx)%ny;
        c2 = 0;
    }
};
