# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 1.33 | 0.0712 | 1.32 | 23 | 6.8e-07 | yes |
| cg-muelu | 5.09 | 3.84 | 5.25 | 20 | 3.8e-07 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 1.41 | 0.148 | 1.4 | 47 | 2.4e-06 | yes |
| cg-muelu | 8.81 | 7.56 | 8.85 | 39 | 1.6e-06 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 0.585 | 0.0542 | 0.663 | 10 | 8.0e-08 | yes |
| cg-muelu | 2.31 | 1.67 | 2.4 | 8 | 9.8e-09 | yes |

