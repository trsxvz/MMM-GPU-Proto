# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.116 | 0.0179 | 0.116 | 15 | 5.4e-06 | yes |
| cg-boomeramg | 0.235 | 0.173 | 0.252 | 30 | 1.2e-07 | yes |
| cg-boomeramg-elasticity | 0.755 | 0.451 | 0.697 | 31 | 1.1e-07 | yes |

## elasticity, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0749 | 0.0408 | 0.0704 | 54 | 6.0e-07 | yes |
| cg-boomeramg-elasticity | 0.0785 | 0.045 | 0.0727 | 53 | 4.9e-07 | yes |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 0.157 | 0.0365 | 0.147 | 25 | 4.4e-05 | yes |
| cg-boomeramg | 0.269 | 0.204 | 0.321 | 47 | 2.5e-07 | yes |
| cg-boomeramg-elasticity | 0.889 | 0.709 | 1.01 | 53 | 4.3e-08 | yes |

## heterogeneous, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.117 | 0.0831 | 0.109 | 94 | 6.6e-07 | yes |
| cg-boomeramg-elasticity | 0.121 | 0.0818 | 0.116 | 94 | 7.4e-07 | yes |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.138 | 0.055 | 0.132 | 13 | 3.7e-08 | yes |
| petsc-cg-gamg | 0.286 | 0.0308 | 0.283 | 8 | 3.7e-07 | yes |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.0541 | 0.014 | 0.0542 | 14 | 4.2e-08 | yes |

