# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.0943 | 0.0912 | 0.0892 | 319 | 9.0e-09 | yes |
| cg-nodal-block-jacobi | 0.102 | 0.101 | 0.104 | 314 | 9.2e-09 | yes |
| cg-boomeramg-elasticity | 0.198 | 0.117 | 0.2 | 30 | 1.2e-07 | yes |
| cg-boomeramg | 0.203 | 0.119 | 0.199 | 30 | 1.2e-07 | yes |
| pcg-boomeramg | 0.205 | 0.136 | 0.215 | 30 | 1.2e-07 | yes |
| ginkgo-cg-jacobi | 0.384 | 0.368 | 0.374 | 623 | 1.0e-08 | yes |
| mkl-cpardiso | 0.774 | 0.0387 | 0.737 | direct | 2.4e-11 | yes |
| petsc-cg-gamg | 1.01 | 0.2 | 0.506 | 177 | 3.5e-06 | yes |
| mumps | 1.22 | 0.0208 | 0.59 | direct | 7.1e-11 | yes |
| strumpack | 1.3 | 0.0434 | 0.813 | 1 | 7.4e-11 | yes |
| petsc-pastix | 1.47 | 0.302 | 1.06 | 1 | 2.8e-11 | yes |
| superlu | 1.73 | 0.0364 | 1.05 | direct | 2.9e-11 | yes |
| gmres-mgr | 4.24 | 3.73 | 3.84 | 362 | 1.3e-05 | yes |
| ginkgo-cg-multigrid | 18.6 | 19.4 | 18.5 | 2000 | 1.3e+02 | no |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-boomeramg | 0.254 | 0.233 | 0.272 | 47 | 2.5e-07 | yes |
| cg-nodal-block-jacobi | 0.257 | 0.259 | 0.26 | 788 | 1.1e-07 | yes |
| cg-boomeramg-elasticity | 0.263 | 0.185 | 0.255 | 47 | 2.5e-07 | yes |
| cg-jacobi | 0.279 | 0.286 | 0.251 | 816 | 1.2e-07 | yes |
| pcg-boomeramg | 0.318 | 0.207 | 0.289 | 47 | 2.5e-07 | yes |
| ginkgo-cg-jacobi | 0.627 | 0.636 | 0.672 | 1043 | 9.9e-09 | yes |
| mkl-cpardiso | 0.745 | 0.0379 | 0.738 | direct | 8.2e-11 | yes |
| petsc-cg-gamg | 1.18 | 0.342 | 0.618 | 284 | 1.3e-06 | yes |
| mumps | 1.18 | 0.0202 | 1.18 | direct | 2.4e-10 | yes |
| strumpack | 1.23 | 0.0432 | 0.8 | 1 | 2.4e-10 | yes |
| petsc-pastix | 1.45 | 0.294 | 1.03 | 1 | 9.9e-11 | yes |
| superlu | 1.71 | 0.0453 | 1.05 | direct | 1.1e-10 | yes |
| gmres-mgr | 11.6 | 11.6 | 12.1 | 1147 | 1.7e-05 | yes |
| ginkgo-cg-multigrid | 20.7 | 19.7 | 19.2 | 2000 | 2.0e+02 | no |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-jacobi | 0.0768 | 0.0716 | 0.0743 | 91 | 7.8e-09 | yes |
| cg-nodal-block-jacobi | 0.0895 | 0.0827 | 0.0852 | 91 | 7.8e-09 | yes |
| pcg-boomeramg | 0.12 | 0.047 | 0.116 | 13 | 3.7e-08 | yes |
| cg-boomeramg | 0.12 | 0.0464 | 0.121 | 13 | 3.7e-08 | yes |
| cg-boomeramg-elasticity | 0.253 | 0.08 | 0.25 | 8 | 4.1e-08 | yes |
| ginkgo-cg-jacobi | 0.393 | 0.36 | 0.409 | 152 | 9.0e-09 | yes |
| ginkgo-cg-multigrid | 0.49 | 0.32 | 0.49 | 38 | 5.9e-09 | yes |
| petsc-cg-gamg | 0.691 | 0.0355 | 0.3 | 8 | 3.7e-07 | yes |
| mkl-cpardiso | 5.63 | 0.206 | 5.62 | direct | 4.7e-14 | yes |
| mumps | 5.66 | 0.0979 | 5.66 | direct | 2.1e-13 | yes |
| strumpack | 5.81 | 0.154 | 5.35 | 1 | 2.1e-13 | yes |
| superlu | 8.23 | 0.293 | 7.72 | direct | 5.6e-14 | yes |
| petsc-pastix | 9.18 | 0.0857 | 15.1 | 1 | 3.1e-13 | yes |

