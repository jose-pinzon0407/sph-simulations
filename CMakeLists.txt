cmake_minimum_required(VERSION 3.11)
include(CMakePrintHelpers)


project(SPH
        VERSION 0.1
        LANGUAGES CXX
        )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


message( STATUS "Project source is: ${CMAKE_CURRENT_SOURCE_DIR}")
set(SPH_CORE "${CMAKE_CURRENT_SOURCE_DIR}/core/")
message( STATUS "Project core is: ${SPH_CORE}")
option(DOUBLE_PRECISION "Either single or double" ON)
option(VECTORIZATION "Use auto vectorization" ON)
option(GPU "Enable CUDA GPU acceleration" OFF)
option(OPENMP "Enable OpenMP" OFF)

set(CONFIG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/config/")
message( STATUS "Project config. directory is: ${CONFIG_DIR}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Settings.h.in 
               ${CONFIG_DIR}/Settings.h
               @ONLY)

if(GPU)
enable_language(CUDA)
set(CMAKE_CUDA_STANDARD 17)
# set(CMAKE_CUDA_ARCHITECTURES native)
add_compile_definitions(GPU=${GPU})
message(STATUS "CUDA acceleration enabled")
endif()

if(OPENMP)
add_compile_definitions(OPENMP=${OPENMP})
message(STATUS "OpenMP Enabled")
find_package(OpenMP)
endif()


set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_compile_options(-O3
                #     -flto
                    )

add_link_options(-flto)

if(DOUBLE_PRECISION)
message(STATUS "Using double precision")
add_compile_definitions(DOUBLE_PRECISION=${DOUBLE_PRECISION})
endif()





if(VECTORIZATION)
message(STATUS "Using auto vectorization")
add_compile_options(;-ftree-vectorize
                     -ffast-math
                     -march=native
                     -fopt-info-vec)
endif()

set(SPH_SOURCES
        ${SPH_CORE}/Kernel.cpp
        ${SPH_CORE}/alg/ForceComputation.cpp
        ${SPH_CORE}/alg/SoAForceComputation.cpp
        ${SPH_CORE}/alg/SortedForceComputation.cpp
        ${SPH_CORE}/compute/ListsForceComputation.cpp
        ${SPH_CORE}/compute/MPListsForceComputation.cpp
        ${SPH_CORE}/alg/NeighborSearch.cpp
        ${SPH_CORE}/container/Container.cpp
        ${SPH_CORE}/container/SoAContainer.cpp
        ${SPH_CORE}/container/Particle.cpp
        ${SPH_CORE}/container/ordering/MortonNeighborSearch.cpp
        ${SPH_CORE}/container/ordering/MortonNeighborLists.cpp
        ${SPH_CORE}/container/ordering/MortonForceComputation.cpp
        ${SPH_CORE}/integrator/Integrator.cpp
        ${SPH_CORE}/utils/Printer.cpp
        ${SPH_CORE}/utils/CheckPointWriter.cpp
        ${SPH_CORE}/utils/Timer.cpp
        ${SPH_CORE}/utils/TimerManager.cpp
        ${SPH_CORE}/utils/InputParser.cpp
        ${SPH_CORE}/alg/NeighborLists.cpp
        ${SPH_CORE}/alg/EncodedListsForceComputation.cpp
)

if(GPU)
        list(APPEND SPH_SOURCES
          ${SPH_CORE}/acceleration/cuContainer.cu
          ${SPH_CORE}/acceleration/cuFindNeighbors.cu
          ${SPH_CORE}/acceleration/cuComputation.cu
          ${SPH_CORE}/acceleration/cuNeighborLists.cu
          ${SPH_CORE}/acceleration/cuComputationWithLists.cu
          )
endif()


add_executable(sim-sph.out
              ${SPH_SOURCES}
              )
if(GPU)
        target_sources(sim-sph.out PRIVATE app/main.cpp)

else()
        target_sources(sim-sph.out PRIVATE app/main.cpp)
endif()



target_include_directories(sim-sph.out PRIVATE ${SPH_CORE})
target_include_directories(sim-sph.out PRIVATE ${CONFIG_DIR})
target_include_directories(sim-sph.out PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/app})

if(OpenMP_CXX_FOUND)
target_link_libraries(sim-sph.out PRIVATE OpenMP::OpenMP_CXX)
endif()

if(GPU)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo --prec-sqrt=false --prec-div=false --use_fast_math")
find_package(CUDAToolkit REQUIRED)
target_link_libraries(sim-sph.out PRIVATE CUDA::cudart)
# set_target_properties(sim-sph.out PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()


message(STATUS "DOUBLE PRECISION:${DOUBLE_PRECISION}")
message(STATUS "VECTORIZATION:${VECTORIZATION}")
message(STATUS "OPENMP:${OPENMP}")
message(STATUS "GPU:${GPU}")

cmake_print_variables(CMAKE_CXX_COMPILER)
cmake_print_variables(CMAKE_BUILD_TYPE)
