#!/bin/bash
# Weak scaling: keep LOCAL_X x LOCAL_Y per MPI rank; grow global grid with ranks.

echo "Weak scaling: multinode with constant local workload"

# --- per-rank workload (fixed) ---
LOCAL_X=4096
LOCAL_Y=4096
N_STEPS=500

# --- hybrid mapping ---
TASKS_PER_NODE=8      # MPI ranks per node
OMP_THREADS=14        # OpenMP threads per rank
CPUS_PER_TASK=${OMP_THREADS}

# --- node list ---
for NODES in 1 2 4 8 16; do
  TOTAL_TASKS=$(( NODES * TASKS_PER_NODE ))

  # ---- choose process grid Px x Py (explicit cases) ----
  # These cases assume TOTAL_TASKS = 8,16,32,64,128 (i.e., TASKS_PER_NODE=8).
  case "${TOTAL_TASKS}" in
    8)    PX=4;  PY=2  ;;   # 1 node  (8 ranks)
    16)   PX=4;  PY=4  ;;   # 2 nodes (16 ranks)
    32)   PX=8;  PY=4  ;;   # 4 nodes (32 ranks)
    64)   PX=8;  PY=8  ;;   # 8 nodes (64 ranks)
    128)  PX=16; PY=8  ;;   # 16 nodes (128 ranks)
    *)
      echo "ERROR: add a case for TOTAL_TASKS=${TOTAL_TASKS} (TASKS_PER_NODE=${TASKS_PER_NODE})" >&2
      exit 1
      ;;
  esac

  # Global sizes so each rank processes LOCAL_X x LOCAL_Y
  GRID_SIZE_X=$(( LOCAL_X * PX ))
  GRID_SIZE_Y=$(( LOCAL_Y * PY ))

  JOB_NAME="weak_${NODES}n_${TOTAL_TASKS}mpi_${OMP_THREADS}omp_${GRID_SIZE_X}x${GRID_SIZE_Y}"

  sbatch \
    --nodes=${NODES} \
    --ntasks=${TOTAL_TASKS} \
    --ntasks-per-node=${TASKS_PER_NODE} \
    --cpus-per-task=${CPUS_PER_TASK} \
    --job-name=${JOB_NAME} \
    --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_THREADS=${OMP_THREADS},JOB_NAME=${JOB_NAME},TOTAL_TASKS=${TOTAL_TASKS},N_TASKS_PER_NODE=${TASKS_PER_NODE} \
    run_scaling.sbatch

  echo "Submitted: nodes=${NODES}, ranks=${TOTAL_TASKS} (grid ${PX}x${PY}), global=${GRID_SIZE_X}x${GRID_SIZE_Y}, local=${LOCAL_X}x${LOCAL_Y}"
done

echo "All Weak Scaling jobs submitted."
