
#include <iostream>

struct grid_1
{
    size_t n[2]{10,10};
    double dx[2]{0,0};
    double start[2]{0,0};
    double end[2]{0,0};

};

struct grid_2
{
   size_t * n;
   double * dx;
   double * start;
   double * end;

};

grid_2 create_grid_2()
{
    grid_2 g;
    g.n = new size_t[2];
    g.dx = new double[2];
    g.start = new double[2];
    g.end = new double[2];

    g.n[0]=10;
    g.n[1]=10;
    g.dx[0]=0;
    g.dx[1]=0;
    g.start[1]=0;    
    g.start[0]=0;
    g.end[1]=0;
    g.end[0]=0;

    
    
    return g;
}

int main( int argc, char ** argv )
{
    grid_1 g;
    std::cout << sizeof(g) << std::endl;
    double * a = new double[10*10];

    #pragma omp target teams distribute  parallel for map(tofrom:a[0:10*10],g.n[0:2],g.dx[0:2],g.start[0:2],g.end[0:2])
    //#pragma omp target teams distribute  parallel for map(tofrom:a[0:10*10],g.n,g.dx,g.start,g.end)
    for(int i=0;i<10;i++)
        for(int j=0;j<10;j++)
        {

            a[ j + i*g.n[0] ] = j + i*g.n[0];
        }
    
     for(int i=0;i<g.n[0];i++)
        for(int j=0;j<g.n[1];j++)
        {

            std::cout << a[ j + i*g.n[0] ]<< std::endl;
        }
}