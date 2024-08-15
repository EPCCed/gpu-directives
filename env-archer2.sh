
module load cpe/23.09
module load PrgEnv-cray
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

run_gpu()
{
    srun --time=0:20:00 --partition=gpu --qos=gpu-shd --gpus=1 $1
}
