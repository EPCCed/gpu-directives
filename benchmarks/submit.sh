#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc

CUR_DIR=$(pwd)
cd ..
source env-archer2.sh 

cd $CUR_DIR

#export CRAY_ACC_DEBUG=3
#export LIBOMPTARGET_KERNEL_TRACE=1

#srun --cpus-per-task=1 --ntasks=1  --hint=nomultithread --distribution=block:block rocprof  --timestamp on   -i input.txt -o results.csv ./build/hip-stream -s 249856

#srun --cpus-per-task=1 --ntasks=1  --hint=nomultithread --distribution=block:block ./build/hip-stream -s 1024 > bandwidth_data.txt

srun --cpus-per-task=1 --ntasks=1  --hint=nomultithread --distribution=block:block python3 run_benchmarks.py  > bandwidth_data.txt