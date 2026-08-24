#pragma once

#include "ConfigVariables.hpp"

void cuKernelUpdatePositions(int blocks, int threads, 
                             int N,
                             Real ubx, Real uby, 
                             Real lbx, Real lby, 
                             Real dt,
                             Real* rx, Real* ry, 
                             Real* vx, Real* vy
                             );

void cuKernelUpdateVelocities(int blocks, int threads, 
                             int N,
                             Real dt,
                             Real* fx, Real* fy, 
                             Real* vx, Real* vy
                             );

void cuCallbackIntegrate(int blocks, int threads,
                         int N,
                         Real ubx, Real uby,
                         Real lbx, Real lby,
                         Real dt,
                         Real* fx, Real* fy,
                         Real* rx, Real* ry,
                         Real* vx, Real* vy);