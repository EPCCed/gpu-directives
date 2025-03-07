
#include "bc.h"

void apply_drichlet_bc(double * field, const grid_t * g, double value)
{

    for(int i=0;i< g->n[0];i++)
        {
            auto k_ghost_left = g->get_index(i,-1);
            auto k_ghost_right = g->get_index(i,g->n[1] );

            field[ k_ghost_left ] = value;
            field[ k_ghost_right ] = value;            
        }
   

   for(int j=0;j<g->n[1];j++)
    {
             auto k_ghost_top = g->get_index(-1,j);
             auto k_ghost_bottom = g->get_index(g->n[0],j );

             field[ k_ghost_bottom ] = value;
             field[ k_ghost_top ] = value; 
    }


}


/**
 * Applies periodic boundary condition by filling ghost cells
*/
void apply_periodic_bc(double * field, const grid_t * g)
{
    //#pragma omp target teams distribute parallel for
    for(int i=0;i< g->n[0];i++)
        {
            auto k_ghost_left = g->get_index(i,-1);
            auto k_right = g->get_index(i,g->n[1] -1 );
            auto k_ghost_right = g->get_index(i,g->n[1] );
            auto k_left = g->get_index(i,0);

            field[ k_ghost_left ] = field[k_right];
            field[ k_ghost_right ] = field[k_right];
            
        }


   //#pragma omp target teams distribute parallel for
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