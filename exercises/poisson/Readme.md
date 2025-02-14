# Poisson solver

This sample code solves the poisson equation $\nabla^2 \phi(x,y) = \rho(x,y)$ on a 2D domain. 
The field is discretized on a regular cartesian grid.
The solution is obtained using a jacobi solver. At every iteration the new value of the field is set according to the values of the old field based on the equation

$\phi^{new}_{i,j}= \frac{1}{2} \left( \frac{\phi( i-1,j) + \phi(i+1,j)}{1+\left(\frac{dx}{dy}\right)^2}\right)  + \frac{1}{2} \left( \frac{\phi( i ,j-1) + \phi(i,j+1)}{1+\left(\frac{dy}{dx}\right)^2}\right) -  \frac{1}{2}\frac{\rho_{i,j}}{\left(\frac{1}{dx}\right)^2 + \left(\frac{1}{dy}\right)^2 }$

We keep iterating until we converge to an approximate solution of the poissson equation.

## Setup your environment for exercises

In the root folder source the correct environment

```bash
cd ../../
source env-archer2.sh
```

## Compiling and running

This folder contains a CPU implementation of the Poisson solver described above. 
Once you have setup your environment, you can compile by typing `make` in this folder.

You can submit a job to the computing nodes using the `submit.sh` slurm batch job script.

```bash
sbatch submit.sh
```

## Exercises

### (1) Offload the main loop to the device using OpenMP offloading

The Jacobi iteration is implemented in the `compute_jacobi` function in `jacobi.cpp`.

Your task is to offload the main loop to the device using OpenMP offloading, making sure to map all the required variables to the device.

You can check the correctness of your solution by using the following command: 

```bash 
python3 check_output.py --compare original_phi0_80000.dat new_phi0_80000.dat
```

* This command should return 0 ( within machine precision). 
* Be sure to save the original output files (phi0_*) for comparison to your completed solution. 
* You will need to be using -X forwarding to use this command. 


### (2) Offload the main loop to the device using OpenMP offloading - hip implementation

The `jacobi_hip.cpp` file contains an alternative implementation of the Jacobi kernel using hip. You cun use this implmentation by compiling with: 

```bash 
make clean; make poisson_rocm
```

This will create a new executable called `poisson_rocm`. Remember to update the name of the executable in the `submit.sh` file.

Your task is to offload this implementation to the device using OpenMP offloading. The hip kernel is submitted by the launch_compute_jacobi_hip and expects device adresses as argument, not host addresses. Use OpenMP directives to pass the correct device adresses to the function. 