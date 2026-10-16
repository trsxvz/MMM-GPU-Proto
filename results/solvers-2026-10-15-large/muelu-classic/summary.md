# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-chebyshev | 0.533 | 0.0253 | 0.563 | 7 | 2.5e-02 | no |
| cg-muelu-jacobi | 0.553 | 0.0399 | 0.583 | 12 | 3.5e-06 | yes |
| cg-muelu-jacobi-reuse | 0.554 | 0.0399 | 0.176 | 12 | 3.5e-06 | yes |
| cg-muelu | 2.41 | 1.86 | 2.42 | 10 | 1.5e-06 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-chebyshev | 0.571 | 0.0668 | 0.605 | 20 | 7.8e-06 | yes |
| cg-muelu-jacobi | 0.59 | 0.0802 | 0.625 | 25 | 6.8e-06 | yes |
| cg-muelu-jacobi-reuse | 0.591 | 0.0802 | 0.216 | 25 | 6.8e-06 | yes |
| cg-muelu | 4.1 | 3.59 | 4.18 | 20 | 3.1e-06 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-chebyshev | 0.313 | 0.067 | 0.302 | 12 | 7.0e-08 | yes |
| cg-muelu-jacobi-reuse | 0.324 | 0.0543 | 0.124 | 10 | 1.3e-07 | yes |
| cg-muelu-jacobi | 0.324 | 0.0538 | 0.317 | 10 | 1.3e-07 | yes |
| cg-muelu | 1.91 | 1.6 | 1.96 | 8 | 1.1e-07 | yes |

