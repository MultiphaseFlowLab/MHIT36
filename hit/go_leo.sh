#!/bin/bash
#SBATCH --account="IscrB_EXCEED"
#SBATCH --job-name="cudec"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test.err
#SBATCH --qos=boost_qos_dbg

#module load nvhpc/25.3
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

chmod 777 binder.sh
mpirun -np 4 --map-by node:PE=8 --rank-by core  ./binder.sh ./mhit36