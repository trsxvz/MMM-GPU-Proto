# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 1.29 | 0.65 | 1.29 | 11 | 1.9e-06 | yes |
| cg-muelu-jacobi-reuse | 1.42 | 0.783 | 1.17 | 15 | 2.8e-06 | yes |
| cg-muelu-jacobi | 1.42 | 0.785 | 1.42 | 15 | 2.8e-06 | yes |
| belos-cg-relaxation | 4.26 | 4.2 | 4.25 | 642 | 1.0e-08 | yes |
| belos-cg-chebyshev | 4.71 | 4.59 | 4.73 | 227 | 9.7e-09 | yes |

## heterogeneous, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 2.13 | 1.49 | 2.12 | 26 | 3.8e-06 | yes |
| cg-muelu-jacobi-reuse | 2.46 | 1.83 | 2.21 | 36 | 3.2e-06 | yes |
| cg-muelu-jacobi | 2.47 | 1.83 | 2.46 | 36 | 3.2e-06 | yes |
| belos-cg-relaxation | 12.2 | 12.2 | 12.2 | 1862 | 1.1e-08 | yes |
| belos-cg-chebyshev | 13.3 | 13.2 | 13.3 | 656 | 9.5e-09 | yes |

## laplace, cpu (2146689 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 1.56 | 0.986 | 1.55 | 9 | 9.7e-08 | yes |
| cg-muelu-jacobi-reuse | 1.76 | 1.19 | 1.39 | 12 | 9.7e-08 | yes |
| cg-muelu-jacobi | 1.76 | 1.19 | 1.75 | 12 | 9.7e-08 | yes |
| belos-cg-relaxation | 3.4 | 3.33 | 3.4 | 183 | 8.3e-09 | yes |
| belos-cg-chebyshev | 3.55 | 3.3 | 3.52 | 63 | 8.6e-09 | yes |

