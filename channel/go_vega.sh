#!/bin/bash
#SBATCH --account="your_account"  # Adjust with your account name
#SBATCH --job-name="cudec"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4

#SBATCH --gres=gpu:4
#SBATCH -p gpu

##SBATCH --gres=gpu:1
##SBATCH -p dev

#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/test_%j.out
#SBATCH --error=slurm_out/test_%j.err

ml purge
ml numactl/2.0.18-GCCcore-13.3.0

ml CMake/3.29.3-GCCcore-13.3.0
ml binutils/2.42-GCCcore-13.3.0

ml NVHPC-HPCX-CUDA12/24.9
ml NVHPC/24.9

CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

# Command used for interactive session --> better testing and debugging:
# salloc -A you_account -p gpu --time=00:20:00 --gpus-per-node=4 --gres=gpu:4 --cpus-per-task=8 --nodes=1 --ntasks-per-node=4 --mem=32GB


chmod 777 binder_vega.sh

mpirun -np 2 --rank-by core  ./binder_vega.sh ./mhit36