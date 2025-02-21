
source ../env-archer2.sh 

set -e
set -x


wget https://github.com/UoB-HPC/BabelStream/archive/refs/tags/v5.0.tar.gz
tar -xvf v5.0.tar.gz 
cp CMakeLists.txt BabelStream-5.0/
rm -r build
mkdir -p build
cd build
CXX=hipcc cmake -DCMAKE_BUILD_TYPE=Debug -DCXX_EXTRA_FLAGS="-x hip -std=c++11 -D__HIP_ROCclr__ --rocm-path=${ROCM_PATH}  -D__HIP_PLATFORM_AMD__ --offload-arch=gfx90a" ../BabelStream-5.0 -DCXX_EXTRA_LINKER_FLAGS=" " -DCXX_EXTRA_LINKER_FLAGS=" "     -DMODEL=hip
make VERBOSE=1