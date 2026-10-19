git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
cmake ..
#enable nvshmem
#cmake -DCUDECOMP_BUILD_EXTRAS=1 -DCUDECOMP_ENABLE_NVSHMEM=1 ..
make -j
