
#--cuda-path=/work/y07/shared/cirrus-software/nvidia/hpcsdk-24.5/Linux_x86_64/24.5/cuda/12.4
module load nvidia/nvhpc-nompi/24.5

CFLAGS=" -g -O3 -fopenmp -fopenmp-targets=nvptx64 -fopenmp-extensions -fopenmp-target-debug  --offload-arch=sm_70 " 
LFLAGS="-L /work/y07/shared/cirrus-software/gcc/10.2.0/lib64 -lstdc++  "
clang++ -c   main.cpp $CFLAGS -o main.o
clang++ main.o $CFLAGS $LFLAGS      -o test