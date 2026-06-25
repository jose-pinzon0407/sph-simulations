#pragma once

#include "math/Math.h"


struct Domain{
    Vec3d l_corner{0,0,0};
    Vec3d h_corner{10,20,0};

    double Volume(){
        double volume=1.0;
        Vec3d lenghts = GetLength();
        //not yet implemented
        std::cout<<"Volume\n";
        return volume;
    }

    Vec3d GetLength(){
        Vec3d lengths;
        for(int i=0;i<lengths.size();++i){
            lengths[i] = h_corner[i]-l_corner[i];
        }

        return lengths;
    }

    bool IsInside(Vec3d r){

        for(int i=0;i<r.size();++i){
            if(!std::isfinite(r[i])) return false;

            if(r[i]<l_corner[i] or r[i]>h_corner[i])
                return false;
            
        }
        return true;
    }
};