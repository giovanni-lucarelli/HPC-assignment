# -------- Build config --------
CC      = mpicc
CFLAGS  = -O3 -Wall -fopenmp -march=$(ARCH)
TARGET  = stencil_mpi

# Architecture target
ARCH    ?= native

# Set ENABLE_OUTPUT=1 to compile with ENABLE_OUTPUT preprocessor flag
ENABLE_OUTPUT ?= 0

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

.PHONY: all clean 

all: $(TARGET)

$(TARGET): $(SRC) $(INC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -o $@ $(SRC) -lm

clean:
	rm -f $(TARGET) *.bin
	rm -rf $(OUT_DIR)
