#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=ompsaxpy
#SBATCH --time=00:01:00
#SBATCH --output=%x-%j.out
#SBATCH --gpus=1
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd

module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

srun -n 1 ./ompsaxpy
