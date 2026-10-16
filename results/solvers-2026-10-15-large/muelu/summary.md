# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity, minimum of the repetitions. Relative residual of the first operator.

## elasticity, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 2.76 | 2.58 | 2.72 | 227 | 9.5e-09 | yes |
| cg-muelu | 5.92 | 2.4 | 5.8 | 20 | 1.2e-06 | yes |
| belos-cg-muelu | 7.26 | 3.79 | 7.16 | 24 | 4.9e-09 | yes |

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 0.826 | 0.317 | 0.839 | 227 | 9.7e-09 | yes |
| cg-muelu | 4.84 | 3.54 | 4.87 | 20 | 4.3e-07 | yes |
| belos-cg-muelu | 6.45 | 5.17 | 6.46 | 23 | 4.5e-09 | yes |

## heterogeneous, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 7.55 | 7.38 | 7.59 | 655 | 9.8e-09 | yes |
| cg-muelu | 8.32 | 4.82 | 8.35 | 41 | 1.6e-06 | yes |
| belos-cg-muelu | 13.8 | 7.73 | 11.1 | 50 | 7.2e-09 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 1.43 | 0.914 | 1.42 | 661 | 9.1e-09 | yes |
| cg-muelu | 7.99 | 6.65 | 8.2 | 38 | 1.8e-06 | yes |
| belos-cg-muelu | 11.5 | 10.2 | 11.5 | 46 | 7.7e-09 | yes |

## laplace, cpu (2146689 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 4.25 | 2.08 | 4.17 | 9 | 1.6e-08 | yes |
| belos-cg-muelu | 5.08 | 2.99 | 5.01 | 10 | 1.3e-09 | yes |

## laplace, gpu (2146689 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| cg-muelu | 2.24 | 1.59 | 2.34 | 8 | 9.8e-09 | yes |
| belos-cg-muelu | 2.88 | 2.3 | 2.98 | 8 | 9.8e-09 | yes |

