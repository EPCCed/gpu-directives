# Makefile for compiling test.f90 with nvfortran

# Compiler
FC = nvfortran

# Compiler flags
FCFLAGS = -cuda -mp=gpu -cudalibs

# Source file
SRC = test.f90

# Target executable
TARGET = a.out

# Rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f $(TARGET) out.log

.PHONY: all clean
