# Poisson solver

This sample code solves the poisson equation $\nabla^2 \phi(x,y) = \rho(x,y)$ on a 2D domain. 
The field is discretized on a regular cartesian grid.
The solution is obtained using a jacobi solver. At every iteration the new value of the field is set according to the values of the old field based on the equation

$\phi^{new}_{i,j}= \frac{1}{2} \left( \frac{\phi( i-1,j) + \phi(i+1,j)}{1+\left(\frac{dx}{dy}\right)^2}\right)  + \frac{1}{2} \left( \frac{\phi( i ,j-1) + \phi(i,j+1)}{1+\left(\frac{dy}{dx}\right)^2}\right) -  \frac{1}{2}\frac{\rho_{i,j}}{\left(\frac{1}{dx}\right)^2 + \left(\frac{1}{dy}\right)^2 }$

We keep iterating until we converge to an approximate solution of the poissson equation.

## Compiling and running

The folder contains a CPU implementation described above. 
First setup your environment

```bash
source ../../env-archer2.sh
```

You can compile by typing `make` in this folder.

You can submit a job to the computing nodes using the `submit.sh` slurm batch job script.

```bash
sbatch submit.sh
```

## Exercise

-   The Jacobi iteration is implemented in the `compute_jacobi` function.
    Offload the main lood to the divice using openmp offloading. Make sure to map all the required variables to the device.
-   Now use a custom mapper to offload the grid object to the device.

