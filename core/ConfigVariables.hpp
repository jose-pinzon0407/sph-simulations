#pragma once

#ifdef DOUBLE_PRECISION
using Real=double;
#else
using Real=float;
#endif

#define H 1
#define D 2
#define TOL 1.0e-12
#define GAS_CONSTANT 30.0
#define RHO_0 25
#define MU 0.01
#define DT 0.01

#include "../config/Settings.h"