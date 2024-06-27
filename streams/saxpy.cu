#include <iostream>

using real_t = double;

// blockDim.x - num threads in a block, .x indicates 1D block labelling
// blockIdx.x - thread index number
// multiplying the above two variables gives start of block
// then add the threadIdx.x offset for the particular thread

__global__ void saxpy_parallel(int n, real_t a, real_t *x, real_t *y)
{


	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i<n)  y[i] = a*x[i] + y[i];


}

void checkCUDAError(const char *msg);

int main()
{
	int N =10;
    int nTrials = 100;
    double tol = 1e-5;


	// allocate vectors on host
	size_t data_size = N * sizeof(real_t);
	real_t* h_x = new real_t[data_size];
	real_t* h_y = new real_t[data_size];

	// allocate device memory
	real_t* d_x; real_t* d_y;

	cudaMalloc( &d_x, data_size);
	cudaMalloc( &d_y, data_size);

    
	for (int i = 0;i<=N-1;i++)
	{
        h_x[i]=1;
        h_y[i]=i;
        
		//std::cout << i << " " <<  h_y[i] << std::endl;
	}

    cudaMemcpy(d_x, h_x, data_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, h_y, data_size, cudaMemcpyHostToDevice);	

	// calculate number of blocks needed for N 
	int nblocks = (N+255)/256;


    std::cout << "Start calculation" << std::endl;

    for(int iTrial=0;iTrial<nTrials;iTrial++)
    {
        saxpy_parallel<<<nblocks,256>>>(N,2.0,d_x,d_y);
    }
    std::cout << "End calculation" << std::endl;
    checkCUDAError("kernel execution calls");


// 	// Copy results back from device memory to host memory
// 	// implicty waits for threads to excute
 	cudaMemcpy(h_y, d_y, data_size, cudaMemcpyDeviceToHost);

	// Check for any CUDA errors
    checkCUDAError("cudaMemcpy calls");

	for (int i = 0;i<=N-1;i++)
	{
        real_t expected = i + nTrials*2;

		if (std::abs( h_y[i] - expected ) > tol )
        {
            std::cout << "Error at " << i << ". Expected "<< expected << " but got " << h_y[i]<<std::endl;
            exit(1);
        }
	}



//   cudaFree(d_x);
//   cudaFree(d_y);

//   delete h_x;
//   delete h_y;

    std::cout << "Completed"<< std::endl;

  return 0;

}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}