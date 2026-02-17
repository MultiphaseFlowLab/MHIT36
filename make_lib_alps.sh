#set environment variables for ALPS build
#uenv 25.7 does not work, but 24.11 does work, so we use that one
#uenv image find prgenv-nvfortran
#uenv image pull prgenv-nvfortran/24.11:v2 
uenv start --view=nvfort prgenv-nvfortran/25.7:v2

git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build

#needed when using 25.7 
export CMAKE_PREFIX_PATH=/user-environment/linux-neoverse_v2/nvhpc-25.7-rcus56e4ntq7bzgxbw6vrdi6celzcpgz/Linux_aarch64/25.7:$CMAKE_PREFIX_PATH

cmake ..
make -j
