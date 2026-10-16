# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| muelu-chebyshev-2 | 0.859 | 0.0246 | 0.835 | 7 | 1.4e+01 | no |
| muelu-jacobi-one-sweep | 0.909 | 0.0613 | 0.883 | 29 | 1.5e-06 | yes |
| muelu-jacobi | 0.913 | 0.0705 | 0.902 | 22 | 1.5e-06 | yes |
| muelu-jacobi-unfiltered | 0.913 | 0.0735 | 0.898 | 23 | 6.4e-07 | yes |
| muelu-jacobi-implicit-transpose | 0.916 | 0.0732 | 0.889 | 23 | 6.0e-07 | yes |
| muelu-jacobi-triple-product | 0.916 | 0.0705 | 0.89 | 22 | 1.6e-06 | yes |
| muelu-jacobi-coarse-10000 | 0.921 | 0.0736 | 0.89 | 23 | 6.8e-07 | yes |
| muelu-chebyshev-2-eigen-30 | 0.922 | 0.0727 | 0.898 | 22 | 9.4e-07 | yes |
| muelu-jacobi-drop-0.02 | 1.07 | 0.145 | 1.01 | 46 | 1.4e-07 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| muelu-chebyshev-2 | 0.966 | 0.133 | 0.942 | 41 | 1.5e-06 | yes |
| muelu-jacobi-one-sweep | 0.97 | 0.123 | 0.943 | 59 | 2.2e-06 | yes |
| muelu-chebyshev-2-eigen-30 | 0.982 | 0.143 | 0.963 | 44 | 2.2e-06 | yes |
| muelu-jacobi | 0.988 | 0.147 | 0.968 | 47 | 1.7e-06 | yes |
| muelu-jacobi-implicit-transpose | 0.988 | 0.147 | 0.969 | 47 | 1.9e-06 | yes |
| muelu-jacobi-coarse-10000 | 0.989 | 0.147 | 0.966 | 47 | 1.8e-06 | yes |
| muelu-jacobi-unfiltered | 0.989 | 0.147 | 0.968 | 47 | 1.9e-06 | yes |
| muelu-jacobi-triple-product | 0.991 | 0.147 | 0.969 | 47 | 1.5e-06 | yes |
| muelu-jacobi-drop-0.02 | 3.6 | 1.75 | 3.57 | 458 | 1.3e-07 | yes |

