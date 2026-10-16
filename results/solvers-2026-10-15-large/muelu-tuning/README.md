# Options of MueLu on the GPU (16/10/2026)

`cases/measure-muelu-options.sh`: variants of the options of MueLu with the
conjugate gradient of MFEM, one process on the GPU (RTX 5080), beam of
421,443 unknowns, homogeneous and with stiff inclusions. mfem-mgis
3a3909cec8 (conversion of the matrices on the GPU), whole machine
reserved, one repetition. Table: `summary.md`. The reference is the
Jacobi smoother of `../muelu-jacobi` (two sweeps damped by 0.6).

- The setup of MueLu (0.82 seconds on the elasticity, 0.06 of which is the
  conversion of the matrix) does not change with the fused triple product
  (`rap: triple product`), the implicit restriction (`transpose: use
  implicit`), the unfiltered matrix of the smoothed prolongator (`sa: use
  filtered matrix`) or a larger coarse problem (`coarse: max size`
  10000). MueLu builds 3 levels, the coarse problem (2000 unknowns at
  most) being solved by KLU2.
- One sweep of the Jacobi smoother: 29 iterations instead of 22, a
  resolution in 0.061 seconds instead of 0.071.
- Dropping the weak connections (`aggregation: drop tol` 0.02) doubles the
  iterations on the elasticity and multiplies them by 10 with the
  inclusions.
- The Chebyshev smoother of degree 2 stops the conjugate gradient on the
  homogeneous beam (not positive definite). With 30 iterations of the
  power method for the largest eigenvalue (10 by default) and a boost
  factor of 1.2, it converges in 22 iterations, as fast as the Jacobi
  smoother.
- A profile by Nsight Systems of the elasticity with the Jacobi smoother
  shows the GPU busy 60 to 120 ms in each interval of 100 ms during the
  setup: the setup is spent in the kernels of MueLu (aggregation,
  products of sparse matrices), not on the host.
