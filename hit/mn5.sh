module purge
module load nvidia-hpc-sdk/24.3
cp Makefile_mn5 Makefile
make clean
make