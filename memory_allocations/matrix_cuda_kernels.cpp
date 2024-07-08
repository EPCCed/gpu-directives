#define BLOCK_SIZE 32

__global__ void matrix_matrix_kernel(double * A, double * B, double * C, int K, int M , int N)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x; // x dimension is the fastest varying
    int i = blockIdx.y * blockDim.y + threadIdx.y;

    if( (i < K) && (j < M))
    {
        //C[i*M + j ]=0;
        
        for(int l=0;l<N;l++)
        {
            C[i*M + j ]+=A[i*N + l]*B[l*M + j];
        }
    }

}


__global__ void matrix_matrix_shared_kernel(double * A, double * B, double * C, int K, int M , int N)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;

    __shared__ double As[BLOCK_SIZE*BLOCK_SIZE];
    __shared__ double Bs[BLOCK_SIZE*BLOCK_SIZE];

    int ai_start =blockIdx.y * blockDim.y;
    int bj_start=blockIdx.x * blockDim.x;

    double tmp=0;


    for( int aj_start=0, bi_start=0;aj_start<M;aj_start+=BLOCK_SIZE,bi_start+=BLOCK_SIZE)
    {

        As[ threadIdx.y*BLOCK_SIZE + threadIdx.x   ] = A[ (ai_start + threadIdx.y)  * N + aj_start + threadIdx.x  ];
        Bs[ threadIdx.y*BLOCK_SIZE + threadIdx.x   ] = B[ (bi_start + threadIdx.y)  * M + bj_start + threadIdx.x  ];
        
        __syncthreads();

        for (int k=0;k<BLOCK_SIZE;k++)
        {
            tmp+= As[threadIdx.y*BLOCK_SIZE  + k  ]*Bs[ k*BLOCK_SIZE + threadIdx.x];
        }

        __syncthreads();

    }

    C[i*M + j ]+=tmp;

}



