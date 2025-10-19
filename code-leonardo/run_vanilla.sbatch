#!/bin/bash
#SBATCH --job-name=stencil_mpi
#SBATCH --output=stencil_mpi_%j.out
#SBATCH --error=stencil_mpi_%j.err

#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=2      # MPI tasks per node
#SBATCH --cpus-per-task=2        # OpenMP threads per MPI rank
#SBATCH --mem=0
#SBATCH --time=00:00:10
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH --exclusive

# --- Load environment modules ---
module purge
module load openmpi/4.1.6--gcc--12.2.0

# --- Build (optional if already built) ---
make clean && make ENABLE_OUTPUT=0

# --- Run ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun ./stencil_mpi -x 1024 -y 1024 -e 3 -E 50 -n 500 -p 0 -o 0
