#!/bin/bash
#SBATCH --account="IscrB_BMO"
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

module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH

chmod 777 binder_leo.sh
mpirun -np 4 --map-by node:PE=8 --rank-by core  ./binder_leo.sh ./mhit36
# many nodes with poisson_dist=1: HCOLL/UCC collectives off (TDMA MPI_Alltoall much faster and more stable)
#mpirun -np 4 --map-by node:PE=8 --rank-by core -x OMPI_MCA_coll_hcoll_enable=0 -x OMPI_MCA_coll_ucc_enable=0 ./binder_leo.sh ./mhit36