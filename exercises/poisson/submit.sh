#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc

WORK_DIR=$(pwd)

source ../../env-archer2.sh 

#export CRAY_ACC_DEBUG=3

cd $WORK_DIR

srun --ntasks=1 --cpus-per-task=1  ./poisson