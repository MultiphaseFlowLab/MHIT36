#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=test.log
#SBATCH --account=IscrB_SONORA
#SBATCH --partition=boost_usr_prod
##SBATCH --qos=boost_qos_dbg
#SBATCH --time=05:30:00

module load gcc
# module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0

srun ./DropCountSize.x
