CXX?= nvc++

CFLAGS= -O2 -cuda  -fmax-errors=1
LFLAGS=  -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -cuda

OBJ = matrix_matrix.o matrix_cuda_kernels.o

all: matrix_matrix


%.o: %.cpp
	$(CXX)  -c -o $@ $< $(CFLAGS)


matrix_matrix: $(OBJ)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

query_gpu: device_query.o
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

#omp_stream: $(OMP_OBJ)
#	$(CXX) $(OMP_FLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean all

clean:
	rm -f *.o *~ core cuda_stream omp_stream