# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity, minimum of the repetitions. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| amesos2-tacho | 1.06 | 0.0168 | 0.965 | direct | 3.9e-11 | yes |
| belos-cg-chebyshev | 1.23 | 1.09 | 1.17 | 113 | 9.1e-09 | yes |
| belos-cg-relaxation | 1.45 | 1.36 | 1.39 | 319 | 9.0e-09 | yes |
| cg-muelu | 3.29 | 0.797 | 3.21 | 23 | 6.4e-07 | yes |
| belos-cg-muelu | 3.62 | 1.11 | 3.52 | 27 | 7.3e-09 | yes |
| belos-cg-frosch | 7.05 | 1.69 | 6.92 | 49 | 8.9e-09 | yes |
| amesos2-klu2 | 37.6 | 0.0701 | 37.6 | direct | 1.2e-10 | yes |

## elasticity, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 0.0885 | 0.0287 | 0.0806 | 114 | 6.7e-09 | yes |
| belos-cg-relaxation | 0.0909 | 0.0322 | 0.0835 | 319 | 9.0e-09 | yes |
| cg-muelu | 0.677 | 0.287 | 0.631 | 22 | 3.8e-07 | yes |
| amesos2-tacho | 0.759 | 0.00302 | 0.715 | direct | 3.4e-11 | yes |
| belos-cg-muelu | 0.843 | 0.447 | 0.791 | 25 | 7.7e-09 | yes |
| belos-cg-frosch | 50 | 0.215 | 48.8 | 2 | 7.0e-11 | yes |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| amesos2-tacho | 1.04 | 0.0167 | 0.962 | direct | 1.3e-10 | yes |
| belos-cg-chebyshev | 2.98 | 2.85 | 2.91 | 300 | 8.3e-09 | yes |
| belos-cg-relaxation | 3.66 | 3.56 | 3.59 | 845 | 9.8e-09 | yes |
| cg-muelu | 4.89 | 2.39 | 4.88 | 71 | 4.6e-07 | yes |
| belos-cg-muelu | 5.5 | 3 | 5.49 | 76 | 7.7e-09 | yes |
| belos-cg-frosch | 7.22 | 1.87 | 7.09 | 54 | 5.7e-09 | yes |
| amesos2-klu2 | 37.9 | 0.0709 | 38 | direct | 5.3e-10 | yes |

## heterogeneous, gpu (56355 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 0.134 | 0.0734 | 0.126 | 299 | 9.3e-09 | yes |
| belos-cg-relaxation | 0.143 | 0.0846 | 0.138 | 845 | 9.7e-09 | yes |
| amesos2-tacho | 0.765 | 0.00309 | 0.706 | direct | 1.2e-10 | yes |
| cg-muelu | 0.904 | 0.497 | 0.844 | 38 | 5.1e-07 | yes |
| belos-cg-muelu | 1.15 | 0.761 | 1.09 | 43 | 9.2e-09 | yes |
| belos-cg-frosch | 50 | 0.218 | 48.9 | 2 | 2.7e-10 | yes |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 0.533 | 0.414 | 0.496 | 32 | 3.4e-09 | yes |
| belos-cg-relaxation | 0.546 | 0.478 | 0.513 | 91 | 7.8e-09 | yes |
| cg-muelu | 1.58 | 0.426 | 1.47 | 8 | 2.0e-08 | yes |
| belos-cg-muelu | 1.72 | 0.577 | 1.61 | 9 | 1.5e-09 | yes |
| amesos2-tacho | 5.94 | 0.0228 | 5.86 | direct | 1.1e-13 | yes |

## laplace, gpu (274625 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.0634 | 0.0139 | 0.0707 | 91 | 7.8e-09 | yes |
| belos-cg-chebyshev | 0.0667 | 0.0129 | 0.0714 | 32 | 3.4e-09 | yes |
| cg-muelu | 0.333 | 0.169 | 0.303 | 7 | 4.8e-08 | yes |
| belos-cg-muelu | 0.415 | 0.265 | 0.393 | 8 | 3.5e-09 | yes |
| amesos2-tacho | 5.4 | 0.00813 | 5.42 | direct | 1.0e-13 | yes |

