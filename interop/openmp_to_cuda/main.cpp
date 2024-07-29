
#include "saxpy_cuda.h"
#include <iostream>

int main(int argc, char * argv)
{
    double * b;
    double * c;
    int n = 100;
    double tol = 1e-9;

    b = new double [n];
    c = new double [n];
    
    for(int i=0;i<n;i++)
    {
        b[i]=i;
        c[i]=0;
    }

    #pragma omp target data map(tofrom:b[0:n],c[0:n])
    {
        #pragma omp target data use_device_addr(c,b)
        {
            saxpy_parallel<<<(n+255)/256,256>>>(n,1.,b,c,1);
        }

    }

    for (int i=0;i<n;i++)
    {
        double expected = b[i];

        if (std::abs(expected - c[i]) > tol )
        {
            std::cout << "Error at " << i << " . Expected " << expected << " instead of " << c[i] << std::endl;
            exit(1);
        }
    }

    std::cout << "Success!" << std::endl;
}