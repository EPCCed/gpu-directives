
#include <iostream>
#include "omp.h"


int main(int argc, char * argv)
{
    double * b;
    double * c;

    int n = 100000;
    int nTrials = 10000;
    bool pinned= false;

    
    double tol = 1e-9;

    b = new double [n];
    //
    if (pinned)
    {
        omp_memspace_handle_t c_memspace = omp_default_mem_space;
        omp_alloctrait_t c_traits[2] = {  { omp_atk_pinned , true   }, {omp_atk_alignment, 128}  } ;
        omp_allocator_handle_t c_alloc = omp_init_allocator(c_memspace,2,c_traits);
        c = (double *) omp_alloc( n * sizeof(double),c_alloc);
        
        //checkCUDAError( cudaHostAlloc( &c, n*sizeof(double)  ,cudaHostAllocDefault));
    }
    else 
    {
        c = new double [n];
    }
    

    for(int i=0;i<n;i++)
    {
        b[i]=i;
        c[i]=0;
    }

    
    #pragma omp target enter data map(to:b[0:n])
    double start=omp_get_wtime();
    for(int iT=0;iT<nTrials;iT++)
    {
        #pragma omp target data map(tofrom:c[0:n])
        {
            #pragma omp target teams distribute parallel for
            for(int i=0;i<n;i++)
            {
                c[i]=1. * c[i] + b[i];
            }
        }

        for(int i=0;i<n;i++)
        {
            c[i]+=1;
        }
    }
    double end=omp_get_wtime();
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