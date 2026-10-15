# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.0369 | 0.0349 | 0.0331 | 319 | 9.0e-09 | yes |
| cg-boomeramg-elasticity | 0.0732 | 0.0414 | 0.0715 | 54 | 6.0e-07 | yes |
| pcg-boomeramg | 0.0758 | 0.0388 | 0.0657 | 54 | 6.0e-07 | yes |
| cg-nodal-block-jacobi | 0.0877 | 0.078 | 0.0728 | 314 | 9.3e-09 | yes |
| cg-boomeramg | 0.111 | 0.0418 | 0.0665 | 54 | 6.0e-07 | yes |
| ginkgo-cg-jacobi | 0.295 | 0.181 | 0.189 | 664 | 9.8e-09 | yes |
| cudss | 0.405 | 0.00416 | 0.298 | direct | 4.3e-11 | yes |
| gmres-mgr | 0.436 | 0.387 | 0.298 | 179 | 1.6e-05 | yes |
| amgx-pcg | 0.447 | 0.34 | 0.438 | 274 | 9.4e-09 | yes |
| cg-amgx | 0.488 | 0.334 | 0.432 | 261 | 5.0e-07 | yes |
| ginkgo-cg-multigrid | 28.1 | 26.9 | 28.5 | 2000 | 1.5e+02 | no |

## heterogeneous, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.0844 | 0.0828 | 0.085 | 816 | 1.2e-07 | yes |
| pcg-boomeramg | 0.107 | 0.0709 | 0.101 | 94 | 6.6e-07 | yes |
| cg-boomeramg-elasticity | 0.111 | 0.0757 | 0.106 | 94 | 6.6e-07 | yes |
| cg-boomeramg | 0.112 | 0.0757 | 0.106 | 94 | 6.6e-07 | yes |
| cg-nodal-block-jacobi | 0.185 | 0.18 | 0.18 | 786 | 1.2e-07 | yes |
| ginkgo-cg-jacobi | 0.35 | 0.293 | 0.305 | 1096 | 9.2e-09 | yes |
| cudss | 0.353 | 0.00416 | 0.303 | direct | 1.4e-10 | yes |
| cg-amgx | 0.587 | 0.487 | 0.563 | 433 | 7.2e-07 | yes |
| amgx-pcg | 0.611 | 0.51 | 0.582 | 457 | 9.2e-09 | yes |
| gmres-mgr | 2.58 | 2.56 | 2.05 | 1192 | 2.0e-05 | yes |
| ginkgo-cg-multigrid | 52.4 | 54.3 | 51.2 | 2000 | 1.8e+02 | no |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.0174 | 0.0157 | 0.0172 | 91 | 7.8e-09 | yes |
| cg-boomeramg | 0.0514 | 0.0148 | 0.0519 | 14 | 4.2e-08 | yes |
| pcg-boomeramg | 0.0546 | 0.0129 | 0.0534 | 14 | 4.2e-08 | yes |
| cg-boomeramg-elasticity | 0.056 | 0.0152 | 0.0554 | 12 | 7.2e-08 | yes |
| cg-nodal-block-jacobi | 0.0654 | 0.0551 | 0.0553 | 91 | 7.8e-09 | yes |
| ginkgo-cg-jacobi | 0.148 | 0.056 | 0.0963 | 91 | 7.8e-09 | yes |
| cg-amgx | 0.157 | 0.0183 | 0.127 | 11 | 1.5e-07 | yes |
| amgx-pcg | 0.157 | 0.0196 | 0.128 | 13 | 6.7e-09 | yes |
| ginkgo-cg-multigrid | 0.444 | 0.35 | 0.441 | 31 | 7.3e-09 | yes |
| cudss | 3.24 | 0.0163 | 3.19 | direct | 1.2e-13 | yes |

