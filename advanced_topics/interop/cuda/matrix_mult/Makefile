# Makefile for compiling test.f90 with nvfortran

# Compiler
FC = nvfortran

# Compiler flags
FCFLAGS = -cuda -mp=gpu -cudalibs -g -Minfo

# Source file
SRC = cuda_matrix_mult.f90

# Target executable
TARGET = a.out

# Rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f $(TARGET) out.log

.PHONY: all clean
