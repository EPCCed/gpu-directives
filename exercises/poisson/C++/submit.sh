#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd


#export CRAY_ACC_DEBUG=3

#srun --ntasks=1 --cpus-per-task=1  ./poisson
srun --ntasks=1 --cpus-per-task=1 rocprof --hip-trace --parallel-kernels ./poisson 