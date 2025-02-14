module load gcc

export CLANG_ROOT=/work/z04/shared/lparisi/clang
export PATH=$CLANG_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CLANG_ROOT/lib:$LD_LIBRARY_PATH
export CPATH=/work/y07/shared/cirrus-software/gcc/10.2.0/include/c++/10.2.0:/work/y07/shared/cirrus-software/gcc/10.2.0/include/c++/10.2.0/x86_64-pc-linux-gnu:$CPATH

