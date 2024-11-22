

/// Contains geometrical and discretisation information

#include <iostream>
#include <cmath>
#include <fstream>

struct grid
{

    inline auto get_index(size_t i, size_t j) const { return (j +1) + (i+1)*(n[1]+2);};

    double x(size_t i) const {return start[0] + dx[0] * i  ;}
    double y(size_t j) const {return start[1] + dx[1] * j  ;}

    size_t size() const { return (n[0]+2) * (n[1] + 2);}
    size_t n[2];
    double dx[2];
    double start[2];
    double end[2];
};


grid make_grid(const double start[2],const double end[2], size_t n[2]  )
{
    struct grid new_grid;
    new_grid.start[0]=start[0];
    new_grid.start[1]=start[1];

    new_grid.end[0]=end[0];
    new_grid.end[1]=end[1];

    new_grid.n[0]=n[0];
    new_grid.n[1]=n[1];

    new_grid.dx[0]=(new_grid.end[0] - new_grid.start[0])/n[0];
    new_grid.dx[1]=(new_grid.end[1] - new_grid.start[1])/n[1];

    return new_grid;
}

/// Makes a step of the jacobi interaction. Compute the new field, given the hold field and the know densitiy field.


void compute_jacobi(double * phi_new, const double * phi_old, const double * rho, const grid * g)
{
    auto nx = g->n[0];
    auto ny = g->n[1];

    auto dx = g->dx[0];
    auto dy = g->dx[1];
    double aspect2 = (dy/dx)*(dy/dx);

    // interior
    for(int i=0;i<nx;i++)
        for( int j=0;j<ny;j++)
        {
            auto k = g->get_index(i,j);

            phi_new[ k ] = 0.5*( (phi_old[k - 1] +  phi_old[k + 1 ])/(1 + 1./aspect2) + (phi_old[k - (ny+2) ] +  phi_old[(k + ny+2) ])/(1 + aspect2) - rho[k] * dx*dx/( 1 + aspect2   ) );
        }
    
};


/// initialise the field with a isotropic gaussian
void init_gaussian(double alpha, double * field, const grid * g)
{
    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            double r2=g->x(i)* g->x(i) + g->y(j)*g->y(j);
            size_t k = g->get_index(i,j);
            field[ k  ] = std::exp( - alpha* r2 );
        }
};

/// Applies periodic boundary condition by filling ghost cells
void apply_periodic_bc(double * field, const grid * g)
{
    for(int i=0;i<g->n[0];i++)
        {
            auto k_ghost_left = g->get_index(i,-1);
            auto k_right = g->get_index(i,g->n[1] -1 );
            auto k_ghost_right = g->get_index(i,g->n[1] );
            auto k_left = g->get_index(i,0);

            field[ k_ghost_left ] = field[k_right];
            field[ k_ghost_right ] = field[k_right];
            
        }
    

   for(int j=0;j<g->n[1];j++)
    {
             auto k_ghost_top = g->get_index(-1,j);
             auto k_ghost_bottom = g->get_index(g->n[0],j );
             auto k_top = g->get_index(0,j);
             auto k_bottom = g->get_index(g->n[0]-1,j );
             

             field[ k_ghost_bottom ] = field[k_top];
             field[ k_ghost_top ] = field[k_bottom];
             
    }
}




/// initialise the field with a constant
void init_constant(double alpha, double * field, const grid * g)
{

    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            size_t k = g->get_index(i,j);
            field[ k  ] = alpha;
        }

};


/// Print field to a binary file
void print_to_file(const double * field, const grid * g,std::string filename)
{
    std::ofstream f;
    size_t nx = g->n[0];
    size_t ny= g->n[1];

    
    f.open(filename,std::ios::binary | std::ios::out);

    f.write( (const char *)&nx, sizeof(size_t)  );
    f.write( (const char *)&ny, sizeof(size_t)  );

    f.write( (const char *)field, sizeof(double)*(g->n[0]+2)*(g->n[1]+2)  );

};

/// initialise the field with the laplacian of a gaussian. Intended for testing
void init_laplacian_gaussian(double alpha, double * field, const grid * g)
{

    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            double r2=g->x(i)* g->x(i) + g->y(j)*g->y(j);
            size_t k = g->get_index(i,j);
            field[ k  ] = -2*( 2 - 2*alpha*r2)* std::exp( - alpha* r2 );
        }

}



/// initialise the field with the laplacian of a gaussian. Intended for testing.



int main(int argc, char ** argv)
{
    size_t nx = 100;
    size_t ny= 100;

    int niterations = 10000;

    double left_box[2]= {-1,-1};
    double right_box[2]= {1,1};
    size_t shape[2] = {nx,ny};

    auto current_grid = make_grid(left_box,right_box,shape);

    auto grid_size = (nx + 2)* (ny + 2);

    double * rho = new double[grid_size];
    double * phi1 = new double[grid_size];
    double * phi2 = new double [grid_size];

    double * phi_old;
    double * phi_new;

    init_laplacian_gaussian(10.0 ,rho,&current_grid);
    init_constant(1.0,phi1,&current_grid);
    init_constant(0.0,phi2,&current_grid);
    apply_periodic_bc(rho, &current_grid);
    apply_periodic_bc(phi1, &current_grid);
    apply_periodic_bc(phi2, &current_grid);

    phi_new = phi1;
    phi_old = phi2;

    for (int i=0;i<niterations;i++)
    {
        std::swap(phi_new,phi_old);
        compute_jacobi(phi_new,phi_old,rho,&current_grid);
        apply_periodic_bc(phi_new, &current_grid);
    }


    print_to_file(phi2,&current_grid,"phi.dat");

}