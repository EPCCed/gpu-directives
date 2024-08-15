CXX?= nvc++

CFLAGS= -O3
LFLAGS=
OMP_FLAGS=  -fopenmp -mp=gpu

CUDA_OBJ = saxpy.o

OMP_OBJ = saxpy_omp.o

all: cuda_stream omp_stream

%.o: %.cu
	$(CXX) -c -o $@ $< $(CFLAGS)

%.o: %.cpp
	$(CXX) $(OMP_FLAGS) -c -o $@ $< $(CFLAGS)

cuda_stream: $(CUDA_OBJ)
	$(CXX) -cuda -o $@ $^ $(LFLAGS) $(LIBS)

omp_stream: $(OMP_OBJ)
	$(CXX) $(OMP_FLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean all

clean:
	rm -f *.o *~ core cuda_stream omp_stream