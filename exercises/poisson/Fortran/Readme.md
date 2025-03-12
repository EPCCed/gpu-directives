# Poisson solver

## Expected performance on MI210 (C++ version )

With 10000 interactions, 500x500 grid for the C++ version

Kernel | Notes | nFields | Time
-- | -- | -- | ---
compute_jacobi | CPU no-opt | 1 |  1.01439 ms/it
compute_jacobi | GPU no-opt + data transfer per iteration| 1 |  1.25582 ms/it
compute_jacobi | GPU no-opt + data transfer once per block | 1 |   0.405076 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) | 1 |   0.405191 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse | 1 |   0.0473415 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse + loop order | 1 |   0.0240578 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + rocm | 1 |   0.0243487 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse + loop order | 10 |   0.247771 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse + loop order + async (1 thread) | 10 |   0.266847 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse + loop order + async(1 thread) | 10 |   0.247771 ms/it
compute_jacobi | GPU no-opt + data transfer once per block (mapper) + collapse + loop order + async (4 threads) | 10 |   0.137915 ms/it


## Expected performance on MI210 ( Fortran version )

With 10000 interactions, 500x500 grid for the Fortran version

Kernel | Notes | nFields | Time
-- | -- | -- | ---
compute_jacobi | CPU no-opt | 1 |  0.48 ms/it
compute_jacobi | GPU no-opt + data transfer per iteration| 1 |  0.8 ms/it
compute_jacobi | GPU no-opt + data transfer once per block | 1 |   0.296 ms/it
compute_jacobi | GPU no-opt + data transfer once per block + collapse | 1 |   0.0458 ms/it
compute_jacobi | GPU no-opt + data transfer once per block + collapse + loop order | 1 |   0.023 ms/it
compute_jacobi | GPU no-opt + data transfer once per block + collapse + loop order | 10 |   0.24 ms/it
compute_jacobi | GPU no-opt + data transfer once per block + collapse + loop order | 10 |   0.6 ms/it

