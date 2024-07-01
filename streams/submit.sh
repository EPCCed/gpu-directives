#!/bin/bash
#
#SBATCH --partition=gpu
#SBATCH --qos=short
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --account=z04
#SBATCH --exclusive

# Load the required modules 
source ../env-cirrus.sh
module load gcc

nsys profile srun ./omp_stream