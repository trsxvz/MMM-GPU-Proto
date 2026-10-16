# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 5.93 | 2.4 | 5.81 | 20 | 1.0e-06 | yes |
| belos-cg-muelu | 7.27 | 3.79 | 7.16 | 24 | 3.4e-09 | yes |

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 4.84 | 3.54 | 4.87 | 20 | 4.8e-07 | yes |
| belos-cg-muelu | 6.49 | 5.23 | 6.5 | 23 | 3.3e-09 | yes |

