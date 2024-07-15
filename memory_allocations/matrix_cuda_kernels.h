#define BLOCK_SIZE 32
#define AK 16
#define BM 16

#define AN 32
#define BN 32


__global__ void matrix_matrix_kernel(double * A, double * B, double * C, int K, int M , int N);
__global__ void matrix_matrix_shared_kernel(double * A, double * B, double * C, int K, int M , int N);
__global__ void matrix_matrix_shared_kernel2(double * A, double * B, double * C, int K, int M , int N);
