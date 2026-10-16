# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.345 | 0.315 | 0.338 | 642 | 9.9e-09 | yes |
| cg-muelu-jacobi-reuse | 0.899 | 0.0736 | 0.203 | 23 | 6.6e-07 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.215 | 0.182 | 0.207 | 183 | 8.3e-09 | yes |
| cg-muelu-jacobi-reuse | 0.305 | 0.0539 | 0.124 | 10 | 8.0e-08 | yes |

