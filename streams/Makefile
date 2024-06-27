CXX?= nvcc

CFLAGS= -O2
LFLAGS=

CUDA_OBJ = saxpy.o

%.o: %.cu
	$(CXX) -c -o $@ $< $(CFLAGS)

cuda_stream: $(CUDA_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~ core cuda_stream