#include "grid.h"


grid_t make_grid(const double start[2],const double end[2], size_t n[2]  )
{
    grid_t new_grid;

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
