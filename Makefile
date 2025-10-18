# -------- Build config --------
CC      = mpicc
CFLAGS  = -O3 -Wall -fopenmp
TARGET  = stencil_mpi

# Set ENABLE_OUTPUT=1 to compile with ENABLE_OUTPUT preprocessor flag
ENABLE_OUTPUT ?= 1

# -------- Paths --------
SRC_DIR     = src
INCLUDE_DIR = include
OUT_DIR     = output_parallel

# -------- Sources --------
SRC = $(SRC_DIR)/stencil_template_parallel.c
INC = $(INCLUDE_DIR)/stencil_template_parallel.h

# If enabled, add the CPP flag
ifeq ($(ENABLE_OUTPUT),1)
  CFLAGS += -DENABLE_OUTPUT
  POST_RUN = mkdir -p $(OUT_DIR); mv -f *.bin $(OUT_DIR)/
else
  POST_RUN = @true
endif

.PHONY: all clean run_serial run_mpi

all: $(TARGET)

$(TARGET): $(SRC) $(INC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -o $@ $(SRC) -lm

# --- Runs ---

# serial run (1 MPI rank)
run_serial: $(TARGET)
	./$(TARGET) -x 256 -y 256 -n 300 -p 1 -o 1
	$(POST_RUN)

# parallel run (MPI + OpenMP)
run_mpi: $(TARGET)
	OMP_NUM_THREADS=3 mpirun -np 4 ./$(TARGET) -x 512 -y 512 -e 3 -E 50 -n 250 -p 1 -o 1
	$(POST_RUN)

clean:
	rm -f $(TARGET) *.bin
	rm -rf $(OUT_DIR)
