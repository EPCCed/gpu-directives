#!/bin/bash
#

#SBATCH --partition=gpu
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00
#SBATCH --account=z04
#SBATCH --exclusive

source ../env-cirrus.sh

export OMP_NUM_THREADS=1
export OMP_DEBUG=true
source ../env-cirrus.sh
source ../clang/env.sh
#srun --ntasks=1 --cpus-per-task=1  ./test

srun --ntasks=1 --cpus-per-task=1 ncu --set full -o omp_sm_report  ./test