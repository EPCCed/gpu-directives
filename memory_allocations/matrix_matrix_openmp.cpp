// Try to compute the matrix product C = A * B
// B : N x M matrix
// A : K x N matrix
// C : K x M matrix
#include <random>
#include <iostream>
#include <mkl_cblas.h>
#include <omp.h>


void matrix_matrix_naive(double * A, double * B, double * C, int K, int M , int N)
{
    #pragma omp target teams distribute parallel for collapse(2)
    for(int i=0;i<K;i++)
        for(int j=0;j<M;j++)
        {
            for(int l=0;l<N;l++)
            {
                C[i*M + j ]+=A[i*N + l]*B[l*M + j];
            }
        }
}

void matrix_matrix_scheduling(double * A, double * B, double * C, int K, int M , int N)
{
    #pragma omp target teams num_teams(K*M/(32*32)) thread_limit(32*32)
    #pragma omp distribute parallel for collapse(2)
    for(int i=0;i<K;i++)
     for(int j=0;j<M;j++)
        {
            for(int l=0;l<N;l++)
            {
                C[i*M + j ]+=A[i*N + l]*B[l*M + j];
            }
        }
}


void set_random_matrix(double * A, size_t N, size_t M,std::mt19937 &rd)
{
    std::uniform_real_distribution<double> dist(0, 1);
    for(size_t i=0;i<N;i++)
        for(size_t j=0;j<M;j++)
        {
            A[i*M + j ]=dist(rd);
        }
}

auto diff_matrix(double * A , double * B,size_t M,size_t N)
{
    double diff=0;
    for(size_t i=0;i<N;i++)
        for(size_t j=0;j<M;j++)
        {
            diff+=std::abs(A[i*M + j ] - B[i*M+j])/(A[i*M + j ]);
        }
    return diff/(N*M);
}

void print_matrix(double * A,size_t N,size_t M)
{
    for(size_t i=0;i< N ;i++)
        {
            for(size_t j=0;j< M ;j++)
                {
                    std::cout << A[i*M + j ] << " ";
                }
            std::cout << std::endl;
        }

}


int main(int argc, char ** argv)
{   

    size_t N=32*4,M=32*4,K=32*8;
    double * A;
    double * B;
    double * C;
    double * C_test;
    const double tol=1e-7;
    int seed = 1039;
    int nTrials=10000;


    A = new double [K * N];
    B = new double [ M * N];
    C = new double [ K *M];
    C_test = new double [ K *M];

    
    std::mt19937 mt( seed );
    set_random_matrix(C,K,M,mt );
    mt = std::mt19937( seed );
    set_random_matrix(C_test,K,M,mt );

    set_random_matrix(A,K,N,mt );
    set_random_matrix(B,N,M,mt );

    for(int i=0;i<nTrials;i++)
    {
        cblas_dgemm(  
            CblasRowMajor,
            CblasNoTrans, CblasNoTrans,
            (int)K,
            (int)M,
            (int)N,
            1e+0,
            A,
            (int)N,
            B,
            (int)M,
            1e+0,
            C_test,
            (int)M
        );
    }

    #pragma omp target enter data map(to:A[0:K*N],B[0:N*M],C[0:K*M],K,M,N)
    auto start = omp_get_wtime();
    for(int i=0;i<nTrials;i++)
    {
        matrix_matrix_scheduling(A, B, C, K, M , N);
    }

    auto end = omp_get_wtime();
    #pragma omp target exit data map(from:A[0:K*N],B[0:N*M],C[0:K*M])

    if(diff_matrix(C,C_test,K,M) > tol )
        {
            std::cout << "Failed! " <<std::endl;
            std::cout << "C " << std::endl;
            print_matrix(C,K,M);
            std::cout << "C ref." << std::endl;
            print_matrix(C_test,K,M);

            exit(1);
        }
    else 
    {
        std::cout << "Success!" << std::endl;
    }

    std::cout << "Time per replica( micros ): " << (end - start)*1e+6/nTrials << std::endl;


}