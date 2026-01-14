#!/bin/bash

echo "Submitting strong scaling jobs (MPI + OpenMP hybrid)"

# Fixed problem size
GRID_SIZE_X=16384
GRID_SIZE_Y=16384
N_STEPS=500

# Fixed OpenMP settings
OMP_THREADS=112
N_TASKS_PER_NODE=1   # MPI ranks per node

# Node list for strong scaling
NODE_LIST="1 2 4 8 16"

for NODES in $NODE_LIST; do
    TOTAL_TASKS=$((NODES * N_TASKS_PER_NODE))
    JOB_NAME="strong_scaling_${NODES}n_${TOTAL_TASKS}mpi_${OMP_THREADS}omp"

    sbatch \
      --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_THREADS=${OMP_THREADS},JOB_NAME=${JOB_NAME},TOTAL_TASKS=${TOTAL_TASKS} \
      --nodes=${NODES} \
      --ntasks-per-node=${N_TASKS_PER_NODE} \
      --cpus-per-task=${OMP_THREADS} \
      --job-name=${JOB_NAME} \
      run_scaling.sbatch
done

echo "All strong scaling jobs submitted."
