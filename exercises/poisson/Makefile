# tool macros
CXX := CC
HIPCXX:= hipcc
CXXFLAGS := -fopenmp -O2 -std=c++17 -D__HIP_PLATFORM_AMD__ 

LFLAGS := -fopenmp -lroctracer64 -lroctx64 -O2
HIPCXXFLAGS:= -O2 -std=c++17 -x hip -D__HIP_ROCclr__ --rocm-path=${ROCM_PATH} -D__HIP_PLATFORM_AMD__ --offload-arch=gfx90a

# compile macros
TARGET_NAME := main
TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

# src files & obj files

OBJ :=	poisson_solver_omp.o

HIP_OBJ := jacobi_kernel.hip.o

# default rule
default: all

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.hip.o: %.cpp
	$(HIPCXX) $(HIPCXXFLAGS) -c -o $@ $<

poisson: $(OBJ) $(HIP_OBJ)
	$(CXX) $(LFLAGS) $(OBJ) $(HIP_OBJ) -o $@

all: poisson

.PHONY: clean

clean:
	rm -f $(OBJ) $(HIP_OBJ)