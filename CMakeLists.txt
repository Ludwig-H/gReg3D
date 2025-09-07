cmake_minimum_required(VERSION 3.18)
project(greg3d LANGUAGES C CXX CUDA)

find_package(CUDAToolkit REQUIRED)

set(CMAKE_BUILD_TYPE Release)

# Enable C++14 for both host and device code
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Build for a modern architecture (sm_70 and newer are supported by recent CUDA)
set(CMAKE_CUDA_ARCHITECTURES 70)
include_directories(
    GDelaunay/Common
    GDelaunay/GDelaunay
    GDelaunay/Main
    GDelaunay/PBA
    )
set(GREG3D_SOURCE_FILES
    GDelaunay/Common/Config.h
    GDelaunay/Common/CudaWrapper.h
    GDelaunay/Common/DtRandom.cpp
    GDelaunay/Common/DtRandom.h
    GDelaunay/Common/Geometry.cu
    GDelaunay/Common/Geometry.h
    GDelaunay/Common/PerfTimer.h
    GDelaunay/Common/predicates.c
    GDelaunay/Common/STLWrapper.h
    GDelaunay/Common/ThrustWrapper.h
    GDelaunay/GDelaunay/GDelaunay.cpp
    GDelaunay/GDelaunay/GDelaunay.h
    GDelaunay/GDelaunay/GDelHost.cu
    GDelaunay/GDelaunay/GDelHost.h
    GDelaunay/GDelaunay/GDelInternal.h
    GDelaunay/GDelaunay/GDelKernels.cu
    GDelaunay/GDelaunay/GDelKernels.h
    GDelaunay/GDelaunay/GDelPredDevice.h
    GDelaunay/GDelaunay/GDelPredKernels.cu
    GDelaunay/GDelaunay/GDelShewchukDevice.h
    GDelaunay/Main/Application.cpp
    GDelaunay/Main/Application.h
    GDelaunay/Main/Main.cpp
    GDelaunay/PBA/pba3D.h
    GDelaunay/PBA/pba3DHost.cu
    GDelaunay/PBA/pba3DKernel.h
    GDelaunay/PBA/Pba.cu
    GDelaunay/PBA/Pba.h
    )
add_executable(greg3d ${GREG3D_SOURCE_FILES})
target_link_libraries(greg3d PRIVATE CUDA::cudart)
