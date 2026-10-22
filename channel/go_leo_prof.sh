#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="cudec_prof"
#SBATCH --time=00:15:00
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH -p boost_usr_prod

# Profiling run: Nsight Systems on a few ranks only (see nsys_rank.sh), e.g. 20-30 time steps (tfin in input.inp)
module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

NP=$((SLURM_NNODES*SLURM_NTASKS_PER_NODE))
# ranks to profile (default in nsys_rank.sh: first, middle and last)
export PROF_RANKS="0 $((NP/2)) $((NP-1))"
mkdir -p prof

chmod 777 binder_leo.sh nsys_rank.sh
mpirun -np ${NP} -x PROF_RANKS --map-by node:PE=8 --rank-by core ./binder_leo.sh ./nsys_rank.sh ./mhit36
