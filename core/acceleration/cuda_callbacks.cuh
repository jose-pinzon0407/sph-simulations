#pragma once
#include "container/SoAContainer.h"
#include "cuParticles.cuh"
#include "Domain.h"

#include<cuda_runtime.h>


/***
 * HIGH ORDER CALLS FROM SIMULATION CONFIG
 * THE TYPE OF ORDER USED IS DECIDED HERE.
 * THIS FILE IS TO BE CALLED FROM THE MAIN.CPP
 * AND SEEN BY THE HOST.
 */




void cuCallbackFindNeighbors(int blocks, int threads, 
                        int N, int* particle_cell_index,
                        int nx, int ny, Real cell_size,
                        int* cell_start, int* cell_end,
                        Real* rx, Real* ry,
                        Real* vx, Real* vy,
                        Real* fx, Real* fy,
                        std::string& order);

void cuCallbackComputeForcesSorted(int blocks, int threads, 
                                 int N, 
                                 Real* rx, Real* ry, 
                                 Real* fx, Real* fy, 
                                 Real* vx, Real* vy, 
                                 int* particle_cell_id, 
                                 int num_cells, int nx, int ny, 
                                 int* start, int* end, 
                                 Real* rho, Real* pressure, 
                                 cuConstants kernel_constants,
                                 std::string& order);

void cuCallbackApproximateDensitySorted(int blocks, int threads, int N, 
                                      Real* rx, Real* ry, 
                                      int* particle_cell_id, 
                                      int num_cells, int nx, int ny, 
                                      int* start, int* end, Real* rho, 
                                      cuConstants kernel_constants,
                                      std::string& order);

void cuCallbackNeighborLists(int blocks, int threads,
                           int N,
                           Real* rx, Real* ry,
                           int* particle_cell_id, 
                           int num_cells, int nx, int ny, 
                           int* cell_start, int* cell_end,
                           int* neighbors_per_particle, int* offset_per_particle,
                           int*& neighbor_list_indices,
                           std::string& order
                          );