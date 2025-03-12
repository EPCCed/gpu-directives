#include "grid.h"


/**
 * Applies Drichlet boundary conditions.
 * @param value Value to be set on all ghost cells
 * @param g Grid definition the space discretization of a field
 * @param field An array of double containing the field values
*/
void apply_drichlet_bc(double * field, const grid_t * g, double value);
