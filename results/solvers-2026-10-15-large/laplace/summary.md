# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## laplace, cpu (2146689 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 1.94 | 1.33 | 1.93 | 14 | 4.3e-08 | yes |
| pcg-boomeramg | 2 | 1.41 | 2 | 14 | 4.3e-08 | yes |
| petsc-cg-gamg | 2.92 | 0.867 | 2.91 | 9 | 3.6e-07 | yes |
| petsc-cg-gamg-true-residual | 3.02 | 0.979 | 3.03 | 11 | 5.2e-09 | yes |
| cg-jacobi | 3.48 | 3.48 | 3.49 | 183 | 8.3e-09 | yes |
| cg-muelu | 4.28 | 2.06 | 4.19 | 9 | 1.6e-08 | yes |
| belos-cg-muelu | 5.1 | 2.96 | 5.02 | 10 | 1.3e-09 | yes |
| ginkgo-cg-multigrid | 6.39 | 3.89 | 6.37 | 52 | 8.6e-09 | yes |
| ginkgo-cg-jacobi | 8.92 | 8.69 | 8.9 | 280 | 9.1e-09 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.19 | 0.184 | 0.183 | 183 | 8.3e-09 | yes |
| pcg-boomeramg | 0.314 | 0.0778 | 0.406 | 16 | 7.4e-08 | yes |
| cg-boomeramg | 0.347 | 0.0885 | 0.415 | 16 | 7.4e-08 | yes |
| cg-amgx-aggregation | 0.485 | 0.24 | 0.552 | 47 | 4.4e-08 | yes |
| ginkgo-cg-jacobi | 0.845 | 0.385 | 0.756 | 183 | 8.3e-09 | yes |
| amgx-pcg | 1 | 0.0926 | 1 | 14 | 4.2e-09 | yes |
| cg-amgx | 1.05 | 0.0864 | 0.99 | 12 | 9.0e-08 | yes |
| ginkgo-cg-multigrid | 1.95 | 1.59 | 1.96 | 46 | 8.0e-09 | yes |

