git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.9/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/25.9/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/25.9/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.9/comm_libs/mpi/man
cmake ..
make -j
