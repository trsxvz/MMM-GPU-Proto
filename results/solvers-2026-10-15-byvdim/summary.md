# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity, minimum of the repetitions. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.116 | 0.0176 | 0.116 | 15 | 5.4e-06 | yes |
| petsc-cg-gamg-true-residual | 0.149 | 0.0233 | 0.119 | 21 | 2.8e-09 | yes |
| cg-boomeramg | 0.203 | 0.124 | 0.203 | 30 | 1.2e-07 | yes |
| cg-boomeramg-elasticity | 0.647 | 0.385 | 0.661 | 31 | 1.1e-07 | yes |

## elasticity, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0743 | 0.0408 | 0.0704 | 54 | 6.0e-07 | yes |
| cg-boomeramg-elasticity | 0.0785 | 0.045 | 0.0702 | 53 | 4.9e-07 | yes |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.131 | 0.0322 | 0.13 | 25 | 4.4e-05 | yes |
| petsc-cg-gamg-true-residual | 0.144 | 0.0437 | 0.134 | 39 | 5.5e-09 | yes |
| cg-boomeramg | 0.263 | 0.204 | 0.276 | 47 | 2.5e-07 | yes |
| cg-boomeramg-elasticity | 0.859 | 0.628 | 0.943 | 53 | 4.3e-08 | yes |

## heterogeneous, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.116 | 0.0814 | 0.109 | 94 | 6.6e-07 | yes |
| cg-boomeramg-elasticity | 0.118 | 0.08 | 0.115 | 94 | 7.4e-07 | yes |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.136 | 0.055 | 0.132 | 13 | 3.7e-08 | yes |
| petsc-cg-gamg | 0.286 | 0.0308 | 0.283 | 8 | 3.7e-07 | yes |
| petsc-cg-gamg-true-residual | 0.312 | 0.0382 | 0.287 | 10 | 4.1e-09 | yes |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0525 | 0.014 | 0.053 | 14 | 4.2e-08 | yes |

