#include <cstdio>
#include <cuda_runtime.h>

void checkCUDAError(cudaError_t err)
{
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString( err) );
        exit(-1);
    }                         
}


__global__ void saxpy_parallel(int n, double a, double *x, double *y,int niters);
