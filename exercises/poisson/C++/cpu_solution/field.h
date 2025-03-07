#ifndef FIELD_H 
#define FIELD_H

#include "grid.h"

struct field_t
{

    void init(grid_t * grid)
    {
        grid=grid;
        data= new double[grid->get_size()]{0};
        size=grid->get_size();
    }
    
    auto get_data() {return data;}
    const auto get_data() const  {return data;}
    size_t get_size(){return size;}


    double * data;
    size_t size;

};

#endif