# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 1.58 | 0.426 | 1.47 | 8 | 2.0e-08 | yes |
| amesos2-tacho | 5.94 | 0.0228 | 5.86 | direct | 1.1e-13 | yes |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 0.333 | 0.169 | 0.303 | 7 | 4.8e-08 | yes |
| amesos2-tacho | 5.4 | 0.00813 | 5.42 | direct | 1.0e-13 | yes |

