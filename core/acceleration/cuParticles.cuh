#pragma once

#include "ConfigVariables.hpp"

struct cuParticles{
    Real* dev_rx=nullptr;
    Real* dev_ry=nullptr;
    Real* dev_rz=nullptr;

    Real* dev_vx=nullptr;
    Real* dev_vy=nullptr;
    Real* dev_vz=nullptr;

    Real* dev_fx=nullptr;
    Real* dev_fy=nullptr;
    Real* dev_fz=nullptr;

    int* cell_start=nullptr;
    int* cell_end=nullptr;
    int* particle_cell_id=nullptr;

    Real* density=nullptr;
    Real* pressure=nullptr;


};
