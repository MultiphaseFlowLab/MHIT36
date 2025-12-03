#!/bin/bash
#SBATCH --account="tra25_openhack"
#SBATCH --job-name="cudec"
#SBATCH --time=00:10:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --qos=boost_qos_dbg
#SBATCH --error=test.err
#### if mapping is on
#SBATCH --cpus-per-task=8


#module load nvhpc/24.3
#module load cuda/12.3
#module load openmpi/4.1.6--nvhpc--24.3
module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/aroccon0/MHIT36_cuDecomp/cuDecomp/build/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/lenzenbe/RE95_256_cuDec/cuDecomp/build/lib:$LD_LIBRARY_PATH
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

# simple run 
#mpirun -n 4 ./mhit36

# profile on single node (all ranks inside one profile)
#nsys profile  -t cuda,nvtx,mpi,openacc  mpirun -np 4 ./mhit36

# profile (one rep for rank)
#mpirun -n 4 nsys profile --trace=cuda,nvtx,mpi -o profile_output_%q{SLURM_PROCID} --stats=true ./mhit36
#mpirun  -n 4 nsys profile -t cuda,nvtx,mpi -o report.$SLURM_LOCALID ./mhit36
#srun -n 4  nsys profile -t cuda,nvtx,mpi --output=nsys_report_rank%t ./mhit36nsys profile --multiprocess=true -t cuda,nvtx,mpi -o report $

# profile + node mapping + nic 
mpirun -np 4 --map-by node:PE=8 --rank-by core nsys profile  -t cuda,nvtx,mpi,openacc  --nic-metrics=true  ./binder.sh ./mhit36

# for nsight compute report
#mpirun -n 4 ncu --kernel-name main_659 --set=full --import-source=yes -o profile -f --launch-skip 3 --launch-count 1 "./mhit36"

# for nsight compute report - all kernels
# mpirun -n 4 ncu --kernel-name regex:main_ --set=full --import-source=yes --launch-skip 70 --launch-count 18 -o reportall.%p ./mhit36

# for nsight compute report - all kernels + mapping + nic
# mpirun -np 4 --map-by node:PE=8 --rank-by core ncu --kernel-name regex:main_ --set=full --import-source=yes --launch-skip 70 --launch-count 18 -o reportall.%p ./binder.sh ./mhit36

