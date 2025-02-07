/**
 * 
 * Initialise the file field with a isotropic gaussian A*exp(-alpha*(x**2 + y**2))
 * @param alpha Width of the gaussian
 * @param field Pointer to field object
 * @param A Gaussian amplitude
 */
void init_gaussian(double alpha, double  * data, grid_t * g,double A=1)
{

    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            double r2=g->x(i)* g->x(i) + g->y(j)*g->y(j);
            size_t k = g->get_index(i,j);
            data[ k  ] = A*std::exp( - alpha* r2 );
        }
};


/**
 * Initialize the field with a constant
 * 
*/
void init_constant(double alpha, double * field, const grid_t * g)
{
    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            size_t k = g->get_index(i,j);
            field[ k  ] = alpha;
        }

};


/**
 * Save the field to a file in binary format
*/
void print_to_file(const double * field, const grid_t * g,std::string filename)
{
    std::ofstream f;
    size_t nx = g->n[0];
    size_t ny= g->n[1];

    
    f.open(filename,std::ios::binary | std::ios::out);

    f.write( (const char *)&nx, sizeof(size_t)  );
    f.write( (const char *)&ny, sizeof(size_t)  );

    f.write( (const char *)field, sizeof(double)*(g->n[0]+2)*(g->n[1]+2)  );

};

/**
 * Initialise the field with the laplacian of a gaussian. Mostly meant for testing.
*/
void init_laplacian_gaussian(double alpha, double * field, const grid_t * g)
{

    for(int i=0;i<g->n[0];i++)
        for( int j=0;j<g->n[1];j++)
        {
            double r2=g->x(i)* g->x(i) + g->y(j)*g->y(j);
            size_t k = g->get_index(i,j);
            field[ k  ] = -2*alpha*( 2 - 2*alpha*r2)* std::exp( - alpha* r2 );
        }

};
