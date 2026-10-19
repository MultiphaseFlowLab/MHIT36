# NVHPC SDK version installed on this machine (change it to match your install)
NVHPC_VERSION=25.9

NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/man
# cuDecomp is built in the repository root (see make_lib_local.sh)
CUDECOMP_DIR=$(cd "$(dirname "$0")/.." && pwd)/cuDecomp/build
export LD_LIBRARY_PATH=$CUDECOMP_DIR/lib:$LD_LIBRARY_PATH
#clean folder output
rm -rf output
mkdir output
cp Makefile_local Makefile
#rm *.dat
make clean
make
mpirun -np 2 ./mhit36 
