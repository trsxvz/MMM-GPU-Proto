# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 1.61 | 0.807 | 1.57 | 19 | 9.6e-06 | yes |
| petsc-cg-gamg-true-residual | 1.77 | 1.05 | 1.76 | 26 | 5.8e-09 | yes |
| cg-boomeramg | 3.92 | 3.35 | 3.82 | 37 | 2.6e-07 | yes |
| pcg-boomeramg | 3.93 | 3.35 | 3.88 | 37 | 2.6e-07 | yes |
| cg-jacobi | 4.65 | 4.67 | 4.63 | 643 | 9.0e-09 | yes |
| cg-nodal-block-jacobi | 4.69 | 4.74 | 4.74 | 635 | 9.4e-09 | yes |
| ginkgo-cg-jacobi | 13.1 | 13.1 | 13.4 | 1353 | 9.3e-09 | yes |
| mumps | 20.1 | 0.2 | 19.8 | direct | 3.7e-10 | yes |
| mkl-cpardiso | 22.2 | 0.588 | 21.8 | direct | 9.8e-11 | yes |

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| pcg-boomeramg | 0.307 | 0.163 | 0.332 | 59 | 9.1e-07 | yes |
| cg-jacobi | 0.336 | 0.33 | 0.334 | 643 | 8.9e-09 | yes |
| cg-boomeramg | 0.346 | 0.188 | 0.357 | 59 | 9.1e-07 | yes |
| cg-nodal-block-jacobi | 0.9 | 0.902 | 0.904 | 635 | 9.4e-09 | yes |
| ginkgo-cg-jacobi | 1.36 | 1.19 | 1.33 | 1347 | 9.8e-09 | yes |
| cg-amgx | 3.85 | 2.53 | 3.84 | 505 | 9.4e-07 | yes |
| amgx-pcg | 4.04 | 2.67 | 4.03 | 539 | 9.4e-09 | yes |
| cudss | 5.84 | 0.0486 | 5.91 | direct | 1.8e-10 | yes |

## heterogeneous, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| petsc-cg-gamg | 2.49 | 1.63 | 2.28 | 40 | 4.5e-05 | yes |
| petsc-cg-gamg-true-residual | 3.03 | 2.35 | 3.03 | 59 | 7.4e-09 | yes |
| cg-boomeramg | 5.12 | 4.66 | 5.09 | 56 | 4.1e-07 | yes |
| pcg-boomeramg | 5.54 | 4.8 | 5.17 | 56 | 4.1e-07 | yes |
| cg-nodal-block-jacobi | 12.6 | 12.6 | 12.6 | 1752 | 1.5e-07 | yes |
| cg-jacobi | 13.6 | 13.3 | 12.5 | 1797 | 1.4e-07 | yes |
| mumps | 19.2 | 0.202 | 19.2 | direct | 1.6e-09 | yes |
| ginkgo-cg-jacobi | 20 | 21.4 | 19.9 | 2000 | 1.2e-01 | no |
| mkl-cpardiso | 20.8 | 0.581 | 21 | direct | 4.7e-10 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| pcg-boomeramg | 0.486 | 0.318 | 0.513 | 109 | 1.5e-06 | yes |
| cg-boomeramg | 0.579 | 0.367 | 0.565 | 109 | 1.5e-06 | yes |
| cg-jacobi | 0.926 | 0.921 | 0.917 | 1797 | 1.4e-07 | yes |
| ginkgo-cg-jacobi | 2.16 | 1.85 | 1.9 | 2000 | 1.7e-02 | no |
| cg-nodal-block-jacobi | 2.52 | 2.48 | 2.47 | 1752 | 1.5e-07 | yes |
| cg-amgx | 5.27 | 4.11 | 5.22 | 898 | 1.6e-06 | yes |
| cudss | 5.4 | 0.0464 | 5.41 | direct | 8.4e-10 | yes |
| amgx-pcg | 5.53 | 4.38 | 5.48 | 956 | 1.0e-08 | yes |

