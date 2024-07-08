CXX?= nvc++

CFLAGS= -O2 -cuda  -fmax-errors=1 -fopenmp
LFLAGS=  -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -cuda -fopenmp -mp=gpu

OBJ = matrix_matrix_cuda.o matrix_cuda_kernels.o
OBJ_OMP = matrix_matrix_openmp.o matrix_cuda_kernels.o

all: matrix_matrix_cuda matrix_matrix_openmp


%.o: %.cpp
	$(CXX)  -c -o $@ $< $(CFLAGS)


matrix_matrix_cuda: $(OBJ)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

matrix_matrix_openmp: $(OBJ_OMP)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)


query_gpu: device_query.o
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

#omp_stream: $(OMP_OBJ)
#	$(CXX) $(OMP_FLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean all

clean:
	rm -f *.o *~ core matrix_matrix_cuda matrix_matrix_openmp