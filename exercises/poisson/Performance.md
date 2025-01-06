# Notes on the performance poisson solver

nFields | shape | hardware | Iterations | Time per it.(ms) | notes
-- | --- | --- | -- | --
1 | 500x500 | 1 CPU | 10 | 0.238957 | base
1 | 500x500 | 1 GPU | 10 | 47.2857 | offload at every iteration
1 | 500x500 | 1 GPU | 10 | 0.296866 | offload only once per block

0.336222

The kernel seems to be limited by latency. It takes about $280 /mu s$ to launch the same empty kernel with openmp.