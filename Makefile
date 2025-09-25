# Compilatore MPI + OpenMP
CC = mpicc
CFLAGS = -O3 -Wall -fopenmp
TARGET = stencil_mpi

# directory
SRC_DIR = src
INCLUDE_DIR = include
OUT_DIR = output_parallel

# sorgenti e header
SRC = $(SRC_DIR)/stencil_template_parallel.c
INC = $(INCLUDE_DIR)/stencil_template_parallel.h

all: $(TARGET)

$(TARGET): $(SRC) $(INC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -o $@ $(SRC) -lm

# --- Esecuzioni ---

# run seriale (1 processo MPI)
run_serial: $(TARGET)
	mkdir -p $(OUT_DIR)
	./$(TARGET) -x 256 -y 256 -n 100 -o 1
	mv -f *.bin $(OUT_DIR)/

# run parallelo (MPI + OpenMP)
run_mpi: $(TARGET)
	mkdir -p $(OUT_DIR)
	export OMP_NUM_THREADS=3; \
	mpirun -np 4 ./$(TARGET) -x 256 -y 256 -e 3 -E 500 -n 500 -p 1 -v 1 -o 1
	mv -f *.bin $(OUT_DIR)/

clean:
	rm -f $(TARGET) *.bin
	rm -rf $(OUT_DIR)
