CC = gcc
CFLAGS = -O3 -Wall
TARGET = stencil

# directory
SRC_DIR = src
INCLUDE_DIR = include
OUT_DIR = output

# sorgenti e header
SRC = $(SRC_DIR)/stencil_template_serial.c
INC = $(INCLUDE_DIR)/stencil_template_serial.h

all: $(TARGET)

$(TARGET): $(SRC) $(INC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -o $@ $(SRC) -lm

run: $(TARGET)
	mkdir -p $(OUT_DIR)
	./$(TARGET) -x 256 -y 256 -n 100 -o 1
	mv plane_*.bin $(OUT_DIR)/

clean:
	rm -f $(TARGET) plane_*.bin
	rm -rf $(OUT_DIR)
