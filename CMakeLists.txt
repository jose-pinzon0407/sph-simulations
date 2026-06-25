cmake_minimum_required(VERSION 3.11)
project(SPH)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message( STATUS "Project source is: ${CMAKE_CURRENT_SOURCE_DIR}")
set(SPH_CORE "${CMAKE_CURRENT_SOURCE_DIR}/core")
message( STATUS "Project core is: ${SPH_CORE}")
option(DOUBLE_PRECISION "Either single or double" ON)

if(DOUBLE_PRECISION)
message(STATUS "Using double precision")
add_compile_definitions(DOUBLE_PRECISION=${DOUBLE_PRECISION})
endif()

add_executable(SPH
              app/main.cpp
              ${SPH_CORE}/Kernel.cpp
              ${SPH_CORE}/container/Container.cpp
              ${SPH_CORE}/container/Particle.cpp
              ${SPH_CORE}/integrator/Integrator.cpp
              )

target_include_directories(SPH PRIVATE ${SPH_CORE})