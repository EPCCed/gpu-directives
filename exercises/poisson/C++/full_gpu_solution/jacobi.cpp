#include "jacobi.h"

/**
 * Makes a step of the jacobi interaction. Compute the new field, given the old field and the known densitiy field.
*/
void compute_jacobi(field_t * phi_new, field_t * phi_old, field_t * rho, int nFields,grid_t * g)
{
    
    #pragma omp parallel
    {
       #pragma omp single
        {
            
            for(int iField=0 ; iField < nFields ; iField++ )
            {

                auto field_phi_new=phi_new[iField].get_data();
                auto field_phi_old=phi_old[iField].get_data();
                auto field_rho=rho[iField].get_data();
                
                #pragma omp task
                    {
                    /* 
                    DATA MAPPING PER STEP IMPLEMENTATION
                    //#pragma omp target data map(tofrom:g->n,g->dx,g->start,g->end,field_phi_new[0:(g->get_size())],field_phi_old[0:(g->get_size())],field_rho[0:(g->get_size())])
                    */
                    
                   #pragma omp target teams distribute parallel for collapse(2)
                   {
                        for(int i=0; i<g->n[0] ; i++)
                            for( int j=0; j<g->n[1] ; j++)
                            {
                                
                                auto k = g->get_index(i,j);

                                auto nx = g->n[0];
                                auto ny = g->n[1];

                                auto dx = g->dx[0];
                                auto dy = g->dx[1];
                                double aspect2 = (dy/dx)*(dy/dx);

                                field_phi_new[ k ] = 0.5*( (field_phi_old[k - 1] +  field_phi_old[k + 1 ])/(1 + 1./aspect2) + (field_phi_old[k - (ny+2) ] +  field_phi_old[(k + ny+2) ])/(1 + aspect2) - field_rho[k] * dx*dx/( 1 + aspect2   ) );
                            }
                    }
                    }
                }    
            
            }

    }

    

};