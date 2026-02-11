#!/bin/bash

git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build

ml purge
ml CMake/3.29.3-GCCcore-13.3.0
ml binutils/2.42-GCCcore-13.3.0

ml NVHPC-HPCX-CUDA12/24.9
ml NVHPC/24.9

cmake ..
make -j$(nproc)
