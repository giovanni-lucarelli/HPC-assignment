#!/bin/bash

echo "Submitting OpenMP thread scaling jobs"

NODES=1
N_TASKS_PER_NODE=1   # one MPI rank
TOTAL_TASKS=1
N_STEPS=500
GRID_SIZE_X=16384
GRID_SIZE_Y=16384

THREAD_LIST="1 2 4 8 16 32 56 84 112"

for OMP_THREADS in $THREAD_LIST; do
    JOB_NAME="omp_scaling_${OMP_THREADS}t"

    sbatch \
      --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_THREADS=${OMP_THREADS},JOB_NAME=${JOB_NAME},TOTAL_TASKS=${TOTAL_TASKS} \
      --nodes=${NODES} \
      --ntasks-per-node=${N_TASKS_PER_NODE} \
      --cpus-per-task=${OMP_THREADS} \
      --job-name=${JOB_NAME} \
      run_threads.sbatch
done

echo "All scaling jobs submitted."
