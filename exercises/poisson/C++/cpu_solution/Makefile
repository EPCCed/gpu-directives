# tool macros
CXX := CC
HIPCXX:= hipcc
CXXFLAGS := -fopenmp -O2 -std=c++17

LFLAGS := -fopenmp -lroctracer64 -lroctx64 -O2
HIPCXXFLAGS:= -O2 -std=c++17 -x hip -D__HIP_ROCclr__ --rocm-path=${ROCM_PATH} -D__HIP_PLATFORM_AMD__ --offload-arch=gfx90a

# compile macros
TARGET_NAME := main
TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

# src files & obj files

ROCM ?= true

OBJ :=	poisson_solver.o timer.o grid.o bc.o

OMP_KERNEL := jacobi.o
HIP_KERNEL := jacobi_kernel.hip.o jacobi_hip.o

# default rule
default: poisson

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.hip.o: %.cpp
	$(HIPCXX) $(HIPCXXFLAGS) -c -o $@ $<

poisson: $(OBJ) $(OMP_KERNEL)
	$(CXX) $(LFLAGS) $(OBJ) $(OMP_KERNEL) -o $@

poisson_rocm: $(OBJ) $(HIP_KERNEL)
	$(CXX) $(LFLAGS) $(OBJ) $(HIP_KERNEL) -o $@

.PHONY: clean

clean:
	rm -f poisson poisson_rocm $(OBJ) $(OMP_KERNEL)  $(HIP_KERNEL) 
