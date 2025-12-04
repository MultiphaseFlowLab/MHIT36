NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.9/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/25.9/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/25.9/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.9/comm_libs/mpi/man
LD_LIBRARY_PATH=/home/milton/MHIT36/cuDecomp/build/lib
#clean folder output
rm -rf output
mkdir output
cp Makefile_local Makefile
#rm *.dat
make clean
make
mpirun -np 2 ./mhit36 
