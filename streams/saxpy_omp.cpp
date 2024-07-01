#include <iostream>
#include "omp.h"
using real_t = double;

int main()
{
	size_t N =256* 30;
    int nTrials = 100;
    int nItPerKernel=1e+5;
    double tol = 1e-5;

	// allocate vectors on host
	real_t* h_x = new real_t[N];
	real_t* h_y = new real_t[N];

    for (size_t i = 0;i<N;i++)
	{
        h_x[i]=1;
        h_y[i]=i;

	}
    

    std::cout << "Start calculation" << std::endl;
    

    #pragma omp target enter data map(to:h_x[0:N],h_y[0:N])
    
    

    for(int iTrial=0;iTrial<nTrials;iTrial++)
    {
        #pragma omp parallel sections
        {
        #pragma omp section
        {
            //std::cout << omp_get_thread_num() << std::endl;
        

        #pragma omp target teams distribute parallel for
              for(int i=0;i<N/2;i++)
              {
                for(int j=0;j<nItPerKernel;j++)
                {
                    h_y[i] = 2*h_x[i] + h_y[i];
                }
              }
        }

        #pragma omp section
        {
            //std::cout << omp_get_thread_num() << std::endl;
            
        #pragma omp target teams distribute parallel for
              for(int i=N/2;i<N;i++)
              {
                for(int j=0;j<nItPerKernel;j++)
                {
                    h_y[i] = 2*h_x[i] + h_y[i];
                }
              }
        }

        }

    }
    #pragma omp target exit data map(from:h_x[0:N],h_y[0:N])
   
    std::cout << "End calculation" << std::endl;


	for (int i = 0;i<N;i++)
	{
        real_t expected = i + nTrials*nItPerKernel*2;
        
		if (std::abs( h_y[i] - expected ) > tol )
        {
            std::cout << "Error at " << i << ". Expected "<< expected << " but got " << h_y[i]<<std::endl;
            
        }
	}


    std::cout << "Completed"<< std::endl;

  return 0;

}

