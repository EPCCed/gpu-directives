module load nvidia/nvhpc-nompi/24.5
module load openmpi/4.1.6-cuda-12.4-nvfortran

export CXX=nvcc
run_gpu()
{
    srun --time=0:20:00 --partition=gpu --qos=short --gres=gpu:1 $1
}