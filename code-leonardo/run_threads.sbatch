#!/bin/bash
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH --time=00:15:00
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

EXEC=./stencil_mpi

module purge
module load openmpi/4.1.6--gcc--12.2.0

# OpenMP settings
export OMP_NUM_THREADS=$OMP_THREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "Running job: $JOB_NAME"
echo "Grid: ${GRID_SIZE_X}x${GRID_SIZE_Y}, Steps: ${N_STEPS}, Threads: ${OMP_THREADS}"

if [[ ${TOTAL_TASKS} -eq 1 ]]; then
    srun -n 1 ${EXEC} -n ${N_STEPS} -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -p 1
else
    srun -n ${TOTAL_TASKS} ${EXEC} -n ${N_STEPS} -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -p 1
fi
