
#include<iostream>
struct A {

    void add_a_gpu() {

        #pragma omp target teams distribute parallel for
        for(int i=0;i<n ; i++)
        {
            c[i]=c[i] + a;
        }
    }

    double a;
    int n;
    double * c; 

};


int main(int argc, char ** argv)
{
    A a;
    double tol = 1e-6;

    a.a = 2;
    a.n = 10;
    a.c = new double [a.n];
    for(int i=0;i<a.n;i++)
    {
        a.c[i]=i;
    }

    #pragma omp declare mapper(A d ) map(d,d.a,d.n,d.c[0:d.n])

    #pragma omp target data map(tofrom:a)
    {
        a.add_a_gpu();
    }


    for(int i=0;i<a.n;i++)
    {
        double expected = ( i + a.a );
        if  (std::abs(a.c[i] - expected) > tol )
        {
            std::cout << "Got " << a.c[i] << "instead of" << expected << std::endl;
            exit(1);
        }
    }


}