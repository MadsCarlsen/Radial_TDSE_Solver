cmake_minimum_required(VERSION 3.31)
project(RadialTDSE)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopt-info-vec-missed=missed.log")
#message(STATUS "Using release flags: ${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(Radial_TDSE_Solver main.cpp
        RadialPropagator.cpp
        RadialPropagator.h
        LinalgRoutines.h
        Types.h
        Fields.h)
