CXX?= nvc++

CFLAGS=   -O2 -cuda  -fmax-errors=1 -mp=gpu -fopenmp -gpu=keep -gpu=lineinfo -gpu=ptxinfo #-gpu=maxregcount:32 -gpu=sm_70
LFLAGS=  -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -gpu=keep -gpu=lineinfo  -cuda -fopenmp -mp=gpu -gpu=ptxinfo -Minfo

OBJ = matrix_matrix_cuda.o matrix_cuda_kernels.o
OBJ_OMP = matrix_matrix_openmp.o matrix_cuda_kernels.o

all: matrix_matrix_cuda matrix_matrix_openmp

OBJ_S = $(OBJ:=.s)
OBJ_OMP_S = $(OBJ_OMP:=.s)


%.o: %.cpp
	$(CXX)  -c -o $@ $< $(CFLAGS)

%.o.s: %.cpp
	$(CXX) -S -o $@ $< $(CFLAGS)

matrix_matrix_cuda: $(OBJ) $(OBJ_S)
	$(CXX) -o $@ $(OBJ) $(LFLAGS) $(LIBS)

matrix_matrix_openmp: $(OBJ_OMP) $(OBJ_OMP_S)
	$(CXX) -o $@ $(OBJ_OMP) $(LFLAGS) $(LIBS)


query_gpu: device_query.o
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

#omp_stream: $(OMP_OBJ)
#	$(CXX) $(OMP_FLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean all

clean:
	rm -f *.o *~ core matrix_matrix_cuda matrix_matrix_openmp