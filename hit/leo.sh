# old modules - ! if you use these, changes in Makefile are needed from -gpu=mem:managed to -gpu=managed
#module load nvhpc/24.3
#module load cuda/12.3
#module load openmpi/4.1.6--nvhpc--24.3
module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
cp Makefile_leonardo Makefile
make clean
make
mkdir -p output
# mpirun -np 4 ./mhit36 
