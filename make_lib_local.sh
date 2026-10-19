# NVHPC SDK version installed on this machine (change it to match your install)
NVHPC_VERSION=25.9

git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/man
cmake ..
make -j
