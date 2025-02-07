#ifndef FIELD_H 
#define FIELD_H

#include "grid.h"

struct field_t
{

    void init(grid_t * grid_)
    {
        grid=grid_;
        data= new double[grid->size()]{0};
        
    }

    auto get_data() {return data;}
    const auto get_data() const  {return data;}
    
    auto get_grid() {return grid;}
    const auto get_grid() const {return grid;}

    double * data;
    grid_t * grid;

};

#endif