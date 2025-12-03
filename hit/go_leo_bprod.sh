#!/bin/bash
#SBATCH --account="IscrB_SONORA"
#SBATCH --job-name="cudec"
#SBATCH --time=00:10:00
#SBATCH --nodes=128    ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8 ## override threads limitation
#SBATCH --output=test.out
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_bprod
#SBATCH --error=test.err

module load nvhpc/24.3
module load cuda/12.3
module load openmpi/4.1.6--nvhpc--24.3
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/aroccon0/MHIT36_cuDecomp/cuDecomp/build/lib:$LD_LIBRARY_$
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/lenzenbe/RE95_256_cuDec/cuDecomp/build/lib:$LD_LIBRARY_P$
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

#export OMP_NUM_THREADS=16
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

mpirun -n 512 ./mhit36
