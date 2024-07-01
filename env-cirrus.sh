module load nvidia/nvhpc-nompi/24.5

export CXX=nvc++


run_gpu()
{
    srun --time=0:20:00 --partition=gpu --qos=short --gres=gpu:1 $1
}