

/// Contains geometrical and discretisation information

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <roctracer/roctx.h>
#include "grid.h"
#include "timer.h"
#include "field.h"
#include "bc.h"
#include "jacobi.h"
#include "tools.h"

int main(int argc, char ** argv)
{
    
    int nFields= 1; // number of equations to solve
    //int niterations = 10000;  // number of time steps
    int niterations = 100000;
    int nIterationsOutput = niterations/5; // Number of iterations between successive outputs

    double left_box[2]= {-1,-1}; // Coordinate of the bottom left corner
    double right_box[2]= {1,1}; // Cooridinat of the top right corner
    size_t shape[2] = { 500 , 500 }; // Grid shape
    
    /**
     * Initialization
    */
   
    std::cout << "Initialise" << std::endl;

    auto current_grid = make_grid(left_box,right_box,shape);
    
    field_t * rho = new field_t[nFields];
    field_t * phi1 = new field_t[nFields];
    field_t * phi2 = new field_t[nFields];

    for (int iField=0 ; iField < nFields; iField++ )
    {
        
        rho[iField].init( &current_grid);
        phi1[iField].init( &current_grid);
        phi2[iField].init(&current_grid);

        init_laplacian_gaussian( 10.0 ,rho[iField].get_data(), &current_grid );
        init_gaussian( 20 ,phi1[iField].get_data(),&current_grid, 0.5 );
        init_gaussian( 20 ,phi2[iField].get_data(),&current_grid, 0.5 );

        apply_drichlet_bc(rho[iField].get_data(), &current_grid, 0);
        apply_drichlet_bc(phi1[iField].get_data(), &current_grid, 0);
        apply_drichlet_bc(phi2[iField].get_data(), &current_grid, 0);

        print_to_file( rho[iField].get_data(), &current_grid, "rho" + std::to_string(iField) + ".dat" );
        print_to_file( phi1[iField].get_data(), &current_grid, "phi" + std::to_string(iField) + "_" + std::to_string(0) + ".dat" );

    }

    field_t * phi_old;
    field_t * phi_new;

    phi_new = phi1;
    phi_old = phi2;

    timer compute_jacobi_timer("compute_jacobi");
    timer total_time_timer("total_time");
    


    total_time_timer.start();

    std::cout << "Start " << niterations << " iterations" << std::endl;
    int i=0;

    while( i<niterations )
    {

        /**
         * Calculations 
        */

        {
            for (int iBlock=0;(iBlock<nIterationsOutput) and (i<niterations);iBlock++)
            {

               std::swap(phi_new,phi_old);
                
                compute_jacobi_timer.start();
               
                compute_jacobi(phi_new,phi_old,rho,nFields);

                compute_jacobi_timer.stop();
            
                i++;

            }
        }
        
            /**
             * Output 
            */
            for(int iField=0;iField<nFields;iField++)
            {
                print_to_file(phi_new[iField].get_data(),&current_grid,"phi" + std::to_string(iField) + "_" + std::to_string(i) + ".dat" );
            }
           
        
    }

    

    total_time_timer.stop();

    std::cout << "Finalize" << std::endl;

    std::cout << total_time_timer << std::endl;
    std::cout << compute_jacobi_timer << std::endl;

    
}
