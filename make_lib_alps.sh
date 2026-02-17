#set environment variables for ALPS build
#uenv 25.7 does not work, but 24.11 does work, so we use that one
#uenv image find prgenv-nvfortran
#uenv image pull prgenv-nvfortran/24.11:v2 
#uenv start prgenv-nvfortran/24.11:v2 --view=nvfort,modules


git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build

cmake ..
make -j
