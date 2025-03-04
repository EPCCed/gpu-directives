#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd

WORK_DIR=$(pwd)
cd ../../
source env-archer2.sh 
cd $WORK_DIR

#export CRAY_ACC_DEBUG=3


srun --ntasks=1 --cpus-per-task=1  ./poisson