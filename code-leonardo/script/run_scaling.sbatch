#!/bin/bash
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH --time=00:05:00
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# nodes / ntasks-per-node / cpus-per-task come dalla submit

EXEC=./stencil_mpi

module purge
module load openmpi/4.1.6--gcc--12.2.0

# --- OpenMP pinning ---
export OMP_NUM_THREADS=${OMP_THREADS:-1}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "Job: $JOB_NAME"
echo "Grid: ${GRID_SIZE_X}x${GRID_SIZE_Y}, Steps: ${N_STEPS}"
echo "Nodes=${SLURM_NNODES}, ranks/node=${SLURM_NTASKS_PER_NODE}, OMP=${OMP_NUM_THREADS}"
echo "TOTAL_TASKS=${TOTAL_TASKS:-$SLURM_NTASKS}"

# ranks per node & total ranks (from sbatch env)
RPN=${SLURM_NTASKS_PER_NODE:-${N_TASKS_PER_NODE:-1}}
NP=${TOTAL_TASKS:-${SLURM_NTASKS:-1}}

# mpirun inside the SLURM allocation, mapped and bound per core
mpirun -np ${NP} \
  --map-by ppr:${RPN}:node:PE=${OMP_NUM_THREADS} \
  --bind-to core \
  ${EXEC} -n ${N_STEPS} -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -p 1 -o 0
