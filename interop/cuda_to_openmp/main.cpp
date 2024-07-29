
#include <iostream>

void checkCUDAError(cudaError_t err)
{
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString( err) );
        exit(-1);
    }                         
}




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

    double * b_device;
    double * c_device;

    checkCUDAError( cudaMalloc( &b_device, n*sizeof(double) ) );
    checkCUDAError(cudaMalloc( &c_device, n*sizeof(double) ) );

    checkCUDAError( cudaMemcpy( b_device, b , n*sizeof(double),cudaMemcpyHostToDevice) );
    checkCUDAError(cudaMemcpy( c_device, c , n*sizeof(double),cudaMemcpyHostToDevice) );    
    
    #pragma omp target teams distribute parallel for simd is_device_ptr(c_device,b_device)
    for(int i=0;i<n;i++)
    {
        c_device[i]=b_device[i];
    }
    
    checkCUDAError(cudaMemcpy( c, c_device , n*sizeof(double),cudaMemcpyDeviceToHost) );
    
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