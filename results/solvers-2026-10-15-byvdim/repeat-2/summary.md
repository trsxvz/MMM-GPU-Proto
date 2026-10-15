# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.121 | 0.0176 | 0.119 | 15 | 5.4e-06 | yes |
| cg-boomeramg | 0.203 | 0.124 | 0.203 | 30 | 1.2e-07 | yes |
| cg-boomeramg-elasticity | 0.647 | 0.385 | 0.661 | 31 | 1.1e-07 | yes |

## elasticity, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0743 | 0.0429 | 0.0716 | 54 | 6.0e-07 | yes |
| cg-boomeramg-elasticity | 0.08 | 0.0459 | 0.0702 | 53 | 4.9e-07 | yes |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.131 | 0.0322 | 0.13 | 25 | 4.4e-05 | yes |
| cg-boomeramg | 0.263 | 0.204 | 0.276 | 47 | 2.5e-07 | yes |
| cg-boomeramg-elasticity | 0.859 | 0.628 | 0.964 | 53 | 4.3e-08 | yes |

## heterogeneous, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.116 | 0.0814 | 0.112 | 94 | 6.6e-07 | yes |
| cg-boomeramg-elasticity | 0.118 | 0.08 | 0.115 | 94 | 7.4e-07 | yes |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.136 | 0.0666 | 0.156 | 13 | 3.7e-08 | yes |
| petsc-cg-gamg | 0.31 | 0.0382 | 0.299 | 8 | 3.7e-07 | yes |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0599 | 0.0165 | 0.0598 | 14 | 4.2e-08 | yes |

