# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 1.33 | 0.0709 | 1.36 | 22 | 1.8e-06 | yes |
| cg-muelu-jacobi-reuse | 1.34 | 0.074 | 0.672 | 23 | 7.6e-07 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 1.41 | 0.149 | 1.41 | 47 | 1.6e-06 | yes |
| cg-muelu-jacobi-reuse | 1.42 | 0.149 | 0.742 | 47 | 2.1e-06 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-jacobi | 0.588 | 0.0541 | 0.676 | 10 | 8.0e-08 | yes |
| cg-muelu-jacobi-reuse | 0.59 | 0.0543 | 0.511 | 10 | 8.0e-08 | yes |

