#!/bin/bash
# Profile only selected MPI ranks with Nsight Systems, the others run normally.
# Usage (inside mpirun, after the binder): ./nsys_rank.sh ./mhit36
# Ranks to profile: PROF_RANKS (space separated), default first, middle and last rank.
NRANKS=${OMPI_COMM_WORLD_SIZE}
RANK=${OMPI_COMM_WORLD_RANK}
PROF_RANKS=${PROF_RANKS:-"0 $((NRANKS/2)) $((NRANKS-1))"}

for r in ${PROF_RANKS}; do
   if [ "${RANK}" -eq "${r}" ]; then
      exec nsys profile -t cuda,nvtx,mpi --mpi-impl=openmpi --cuda-memory-usage=false \
           -f true -o prof/prof_${NRANKS}_rank${RANK} "$@"
   fi
done
exec "$@"
