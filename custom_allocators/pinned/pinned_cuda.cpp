
#include <iostream>
#include "saxpy_cuda.h"
#include "omp.h"
#include <stdexcept>
#include <cuda_runtime.h>
int main(int argc, char * argv)
{
    double * b;
    double * c;

    int n = 100000;
    int nTrials = 10000;
    std::string pinning="omp";


    double tol = 1e-9;

    b = new double [n];
    //
    if (pinning == "omp")
    {
        //
        omp_memspace_handle_t c_memspace = omp_default_mem_space;
        omp_alloctrait_t c_traits[2] = {  { omp_atk_pinned , true   }, {omp_atk_alignment, 128}  } ;
        omp_allocator_handle_t c_alloc = omp_init_allocator(c_memspace,2,c_traits);
        c = (double *) omp_alloc( n * sizeof(double),c_alloc);
    }
    else if (pinning == "cuda")
    {
        checkCUDAError( cudaHostAlloc( &c, n*sizeof(double)  ,cudaHostAllocDefault));
    }
    else if (pinning == "none")
    {
        c = new double [n];
    }
    else 
    {
        throw std::invalid_argument("Unkown pinning");

    }
    

    for(int i=0;i<n;i++)
    {
        b[i]=i;
        c[i]=0;
    }

    double * b_device;
    double * c_device;
    checkCUDAError( cudaMalloc( &b_device, n*sizeof(double) ) );
    checkCUDAError(cudaMalloc( &c_device, n*sizeof(double) ) );


    checkCUDAError( cudaMemcpy( b_device, b , n*sizeof(double),cudaMemcpyHostToDevice) );

    double start=0,end=0;
    start=omp_get_wtime();

    for(int iT=0;iT<nTrials;iT++)
    {
        checkCUDAError(cudaMemcpy( c_device, c , n*sizeof(double),cudaMemcpyHostToDevice) );    

        saxpy_parallel<<<(n+256)/256,256>>>(n, 1, b_device, c_device, 1 );
        cudaDeviceSynchronize();

        checkCUDAError( cudaGetLastError() );
        checkCUDAError( cudaMemcpy( c, c_device , n*sizeof(double),cudaMemcpyDeviceToHost) );
    

        for(int i=0;i<n;i++)
        {
            c[i]+=1;
        }
    }
    end=omp_get_wtime();
    std::cout << end - start << std::endl;

    
    for (int i=0;i<n;i++)
    {
        double expected = nTrials*(b[i] + 1);

        if (std::abs(expected - c[i]) > tol )
        {
            std::cout << "Error at " << i << " . Expected " << expected << " instead of " << c[i] << std::endl;
            exit(1);
        }
    }
    std::cout << "Success!" << std::endl;
}