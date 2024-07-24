
#include<iostream>
#include <omp.h> 

#define BLOCK_SIZE 32

int main(int argc, char ** argv)
{
    omp_allocator_handle_t my_allocator;
    omp_alloctrait_t traits[1];

    traits[0].key=omp_atk_access;
    traits[0].value=omp_atv_cgroup;
    omp_memspace_handle_t memspace { omp_default_mem_space };

    my_allocator = omp_init_allocator( memspace , 1, traits );
    int n = BLOCK_SIZE;

    //double * c = (double *)omp_alloc( n * sizeof(double),my_allocator);
    double sum=0;
    double c[BLOCK_SIZE];
    
    #pragma omp target data  map(tofrom:sum)
    {
        
        #pragma omp target teams num_teams(1) reduction(+:sum) shared(my_allocator) private(c)   uses_allocators(omp_pteam_mem_alloc) allocate(omp_pteam_mem_alloc:c)
        {
            //double c[BLOCK_SIZE];
            //#pragma omp allocate(c) allocator(omp_pteam_mem_alloc)
            //#pragma omp allocate(c) allocator(my_allocator)

            #pragma distribute parallel for 
            for(int i=0;i<BLOCK_SIZE;i++)
            {
                c[i]=i;
            }
            
            for(int i=0;i<n;i++)
            {
                sum+=c[i];
            }

        }
    

    }

    std::cout << sum << std::endl;
    

}