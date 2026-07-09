cmake_minimum_required(VERSION 3.11)
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

set(default_build_type "Release")
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

add_executable(SPH_SoA
              app/main.cpp
              ${SPH_CORE}/Kernel.cpp
              ${SPH_CORE}/alg/ForceComputation.cpp
              ${SPH_CORE}/alg/SoAForceComputation.cpp
              ${SPH_CORE}/container/Container.cpp
              ${SPH_CORE}/container/SoAContainer.cpp
              ${SPH_CORE}/container/Particle.cpp
              ${SPH_CORE}/integrator/Integrator.cpp
              ${SPH_CORE}/utils/Printer.cpp
              ${SPH_CORE}/utils/CheckPointWriter.cpp
              ${SPH_CORE}/utils/Timer.cpp
              )

target_include_directories(SPH_SoA PRIVATE ${SPH_CORE})