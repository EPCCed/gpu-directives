#!/bin/bash

#SBATCH --job-name=poisson
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd


WORK_DIR=$(pwd)
cd ../../
source env-archer2.sh
cd $WORK_DIR


module use /mnt/lustre/a2fs-work4/work/z19/shared/rocm/install/3.0.0/modulefiles
module load rocprofiler-compute

#srun  --ntasks=1 --cpus-per-task=1 ./poisson
EXE=poisson

rm -rf workloads

srun  --ntasks=1 --cpus-per-task=1 rocprof-compute profile -n $EXE -- ./$EXE

rocprof-compute analyze -p workloads/poisson/MI200/ > report.txt
