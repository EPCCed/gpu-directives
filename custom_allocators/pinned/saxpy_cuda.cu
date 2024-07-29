
__global__ void saxpy_parallel(int n, double a, double *x, double *y,int niters)
{

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i<n)  
    {
        for(int j=0;j<niters;j++)
        {
            y[i] = a*x[i] + y[i];
        }
    }    

}
