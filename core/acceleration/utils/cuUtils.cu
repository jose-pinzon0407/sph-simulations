

#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>

#include "cuUtilsLauncher.cuh"


struct MaximumVelocity{
    __device__
    Real operator()(const thrust::tuple<Real,Real>& tp)const{
        Real vx = thrust::get<0>(tp);
        Real vy = thrust::get<1>(tp);

        return (vx*vx+vy*vy);
    }

};

struct KineticEnergy{
    __device__
    Real operator()(const thrust::tuple<Real,Real>& tp)const{
        Real vx = thrust::get<0>(tp);
        Real vy = thrust::get<1>(tp);

        return 0.5*(vx*vx+vy*vy);
    }

};

void cuKernelMeasurementsPerStep(int blocks, int threads, 
                                int N,
                                Real& max_vel, Real& kinetic_energy, 
                                Real* rx, Real* ry, 
                                Real* vx, Real* vy){

    auto begin = thrust::make_zip_iterator(thrust::make_tuple(vx,vy));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(vx,vy))+N;

    max_vel = thrust::transform_reduce(
        thrust::device,
        begin,
        end,
        MaximumVelocity{},
        Real{0},
        thrust::maximum<Real>()
    );


    kinetic_energy = thrust::transform_reduce(
        thrust::device,
        begin,
        end,
        KineticEnergy{},
        Real{0},
        thrust::plus<Real>()
    );

}