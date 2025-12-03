#!/bin/bash
#SBATCH --ntasks=16  # --ntasks=8 when used two nodes
#SBATCH --account="ehpc244"
#SBATCH --job-name="cudec"
#SBATCH --time=00:5:00
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=20
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --qos=acc_debug
#SBATCH --error=test.err

module load nvidia-hpc-sdk/24.3
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH


mpirun -n 8 ./mhit36