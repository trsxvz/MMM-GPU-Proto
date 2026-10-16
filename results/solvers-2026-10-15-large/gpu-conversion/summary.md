# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.385 | 0.314 | 0.38 | 642 | 9.9e-09 | yes |
| cg-muelu-chebyshev | 0.861 | 0.0215 | 0.84 | 6 | 1.2e+01 | no |
| cg-muelu-jacobi | 0.91 | 0.0704 | 0.896 | 22 | 1.9e-06 | yes |
| cg-muelu-jacobi-reuse | 0.914 | 0.0705 | 0.24 | 23 | 8.7e-07 | yes |
| cg-muelu-jacobi-unsmoothed | 0.944 | 0.185 | 0.919 | 65 | 4.3e-07 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu-chebyshev | 0.968 | 0.13 | 0.941 | 40 | 2.0e-06 | yes |
| belos-cg-relaxation | 0.983 | 0.909 | 0.977 | 1862 | 1.0e-08 | yes |
| cg-muelu-jacobi | 0.99 | 0.147 | 0.971 | 47 | 1.7e-06 | yes |
| cg-muelu-jacobi-reuse | 0.993 | 0.147 | 0.318 | 47 | 1.8e-06 | yes |
| cg-muelu-jacobi-unsmoothed | 1.09 | 0.336 | 1.07 | 119 | 7.9e-07 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.218 | 0.182 | 0.212 | 183 | 8.3e-09 | yes |
| cg-muelu-chebyshev | 0.289 | 0.0621 | 0.272 | 11 | 2.1e-07 | yes |
| cg-muelu-jacobi | 0.305 | 0.0538 | 0.29 | 10 | 8.0e-08 | yes |
| cg-muelu-jacobi-reuse | 0.305 | 0.0538 | 0.128 | 10 | 8.0e-08 | yes |
| cg-muelu-jacobi-unsmoothed | 0.373 | 0.193 | 0.351 | 40 | 4.5e-08 | yes |

