# Linear solvers

Times in seconds: setup (`SetOperator`) plus first resolution, following resolutions (mean), setup plus resolution for a second operator of the same sparsity. Relative residual of the first operator.

## elasticity, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.0901 | 0.071 | 0.0818 | 319 | 9.0e-09 | yes |
| belos-cg-chebyshev | 0.0943 | 0.0723 | 0.0853 | 114 | 6.7e-09 | yes |
| cg-muelu | 0.14 | 0.0367 | 0.127 | 11 | 1.2e-06 | yes |
| belos-cg-muelu | 0.141 | 0.0343 | 0.124 | 14 | 8.9e-09 | yes |
| cg-muelu-jacobi-reuse | 0.144 | 0.0389 | 0.106 | 14 | 1.3e-06 | yes |
| cg-muelu-jacobi | 0.144 | 0.0402 | 0.156 | 14 | 1.3e-06 | yes |
| amesos2-tacho | 1.04 | 0.0175 | 0.956 | direct | 3.9e-11 | yes |
| belos-cg-frosch-tacho | 4.95 | 2.16 | 4.79 | 49 | 8.9e-09 | yes |
| belos-cg-frosch | 7.19 | 1.69 | 6.88 | 49 | 8.9e-09 | yes |

## heterogeneous, cpu (56355 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-muelu | 0.168 | 0.0634 | 0.155 | 27 | 5.4e-09 | yes |
| cg-muelu | 0.173 | 0.0669 | 0.16 | 21 | 1.8e-06 | yes |
| cg-muelu-jacobi | 0.181 | 0.0774 | 0.169 | 29 | 1.4e-06 | yes |
| cg-muelu-jacobi-reuse | 0.183 | 0.0766 | 0.136 | 29 | 1.4e-06 | yes |
| belos-cg-relaxation | 0.204 | 0.191 | 0.199 | 845 | 9.8e-09 | yes |
| belos-cg-chebyshev | 0.209 | 0.191 | 0.202 | 300 | 8.4e-09 | yes |
| amesos2-tacho | 1.01 | 0.0165 | 0.934 | direct | 1.3e-10 | yes |
| belos-cg-frosch-tacho | 5.22 | 2.4 | 4.95 | 54 | 5.7e-09 | yes |
| belos-cg-frosch | 7.25 | 1.87 | 7.25 | 54 | 5.7e-09 | yes |

## laplace, cpu (274625 unknowns, 8 processes)

| configuration | setup + first | solve | second operator | iterations | residual | converged |
|---|---|---|---|---|---|---|
| belos-cg-relaxation | 0.067 | 0.047 | 0.0568 | 91 | 7.8e-09 | yes |
| belos-cg-chebyshev | 0.0783 | 0.0532 | 0.0676 | 32 | 3.5e-09 | yes |
| belos-cg-muelu | 0.12 | 0.045 | 0.108 | 9 | 4.3e-09 | yes |
| cg-muelu | 0.129 | 0.0516 | 0.117 | 8 | 4.7e-08 | yes |
| cg-muelu-jacobi | 0.152 | 0.072 | 0.139 | 12 | 2.5e-08 | yes |
| cg-muelu-jacobi-reuse | 0.153 | 0.0719 | 0.111 | 12 | 2.5e-08 | yes |
| amesos2-tacho | 5.9 | 0.0222 | 5.83 | direct | 1.1e-13 | yes |

