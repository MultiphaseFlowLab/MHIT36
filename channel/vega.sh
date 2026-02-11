#!/bin/bash

ml purge
ml CMake/3.29.3-GCCcore-13.3.0
ml binutils/2.42-GCCcore-13.3.0


ml NVHPC-HPCX-CUDA12/24.9
ml NVHPC/24.9

cp Makefile_vega Makefile
make clean
make -j$(nproc)
mkdir -p output
