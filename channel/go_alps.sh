#!/bin/bash
#SBATCH --account=lp146
#SBATCH --job-name=example-%j
#SBATCH --time=00:5:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4               # total ranks
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --ntasks=8

CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

srun  mhit36
