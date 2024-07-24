git clone https://github.com/llvm/llvm-project.git
cd llvm-project
git checkout tags/llvmorg-18.1.8
cd ../
mkdir build
cd build
module load gcc
module load cmake
export CPATH=/work/y07/shared/cirrus-software/gcc/10.2.0/include/c++/10.2.0:/work/y07/shared/cirrus-software/gcc/10.2.0/include/c++/10.2.0/x86_64-pc-linux-gnu:$CPATH
cmake ../llvm-project/llvm -DLLVM_ENABLE_PROJECTS="clang" -DLLVM_ENABLE_RUNTIMES="openmp"  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/work/z04/shared/lparisi/clang

