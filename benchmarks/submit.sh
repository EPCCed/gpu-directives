#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc

source ../env-archer2.sh 
#export CRAY_ACC_DEBUG=3
#export LIBOMPTARGET_KERNEL_TRACE=1

#srun --cpus-per-task=1 --ntasks=1  --hint=nomultithread --distribution=block:block rocprof  --timestamp on   -i input.txt -o results.csv ./build/hip-stream -s 249856

srun --cpus-per-task=1 --ntasks=1  --hint=nomultithread --distribution=block:block ./build/hip-stream -s 134217728