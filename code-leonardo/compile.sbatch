#!/bin/bash

#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH --time=00:05:00

#SBATCH --job-name=compile_stencil
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0

module purge
module load openmpi/4.1.6--gcc--12.2.0

make clean && make ENABLE_OUTPUT=0
