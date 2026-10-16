# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 2.76 | 2.58 | 2.72 | 227 | 9.5e-09 | yes |

## elasticity, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 0.826 | 0.317 | 0.839 | 227 | 9.7e-09 | yes |

## heterogeneous, cpu (421443 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 7.55 | 7.38 | 7.59 | 655 | 9.8e-09 | yes |

## heterogeneous, gpu (421443 unknowns, 1 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-chebyshev | 1.43 | 0.914 | 1.42 | 661 | 9.1e-09 | yes |

