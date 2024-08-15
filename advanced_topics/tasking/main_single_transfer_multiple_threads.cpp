#include <iostream>

int main(int argc, char ** argv)
{
    int n=10;
    int m = 100000;
    double alpha=10;
    int nTrials=10000;
    double tol = 1e-7;
    

    double ** xs;
    double ** ys;

    xs = new double* [n];
    ys = new double* [n];

    for (int i=0;i<n;i++)
    {
        xs[i] = new double [m];
        ys[i] = new double [m];
        double * xsi = xs[i];
        double * ysi = ys[i];
        #pragma omp target enter data map(to:xsi[0:m],ysi[0:m])

    }

    std::cout << "Start. " << std::endl;
    
    #pragma omp parallel for
            for ( int islice = 0 ; islice < n ; islice++)
                {
                double * xsi = xs[islice];
                double * ysi = ys[islice];

                
                #pragma omp target teams distribute parallel  for nowait depend(inout:xsi[0:m],ysi[0:m])
                for(int j=0;j<m;j++)
                {
                    xsi[  j  ]=j;
                    ysi[ j ]=0;
                }
                
                #pragma omp target teams distribute  parallel for nowait depend(inout:xsi[0:m],ysi[0:m])
                for(int j=0;j<m;j++)
                {
                    for( int iT=0;iT<nTrials;iT++)
                    {
                        ysi[j]=ysi[j ] + alpha*xsi[ j];
                    }
                    
                }

                }
    
    #pragma omp taskwait

    for (int islice = 0 ; islice < n ; islice++)
    {
        double * xsi = xs[islice];
        double * ysi = ys[islice];

        #pragma omp target exit data map(from:xsi[0:m],ysi[0:m])



    }

    double sum=0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            sum+=ys[i][j];
        }
        
    }

    double expected = m*1.*(m-1)/2 * n * nTrials*alpha ;

    if (std::abs(expected - sum ) > tol )
    {
        std::cout << "Error: sum " << sum <<  " differs from "<< expected <<  std::endl;
        exit(1);
    }
    
    delete xs;
    delete ys;

    std::cout << "Success." << std::endl;


}