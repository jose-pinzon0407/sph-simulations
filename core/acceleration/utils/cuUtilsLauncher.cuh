#pragma once

#include "ConfigVariables.hpp"

void cuKernelMeasurementsPerStep(int blocks, int threads, 
                                 int N,
                                 Real& mv, Real& ke, 
                                 Real* rx, Real* ry, 
                                 Real* vx, Real* vy);