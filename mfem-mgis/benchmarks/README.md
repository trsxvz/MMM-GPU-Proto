# Benchmarks section


Please, refer to the docuementation website for more information: https://thelfer.github.io/mfem-mgis/developer_guide/benchmark.html

It includes:

- Description of the profiling tools
- Benchmarks:
  - RVE of a MOX simulation from example 7 of `mfem-mgis-examples`: https://github.com/latug0/mfem-mgis-examples/tree/master/ex7
  - Overhead between MFEM and MFEM-MGIS 

## Overhead between MFEM and MFEM-MGIS

This benchmark aims to evaluate the time overhead between `MFEM` and `MFEM/MGIS` for a simple test case involving an elastic behavior law.

Details:

```
    Solver: HyprePCG

    Preconditionner: HypreDiagScale (Jacobi)

    Tolerance: 1e-14

    Boundary conditions: we impose an uniform dirichlet condition U = (0,0,0) at left and U = (0,0,1) at right.
```

Note that due to the small size of the input mesh, domain decomposition is limited to 16 subdomains (MFEM-MGIS only performs parallel refinement). To improve this benchmark, consider using sequential refinement or a finer mesh.

Run this example on Topaze supercomputer:

```
ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5
ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMMGISLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5
```

## Linear solvers

`LinearSolversBenchmark` measures a linear solver of the factory of linear
solvers (option `--solver`, its parameters being given by
`--solver-parameters` in a text syntax described in the source) on an
assembled system: the Laplace equation on a cube (`--problem laplace`) or
the linear elasticity of a cantilever beam (`--problem elasticity`), with
`--size` elements across the section, cubic inclusions whose coefficients
are multiplied by `--contrast`, and the unknowns ordered by nodes or by
vector dimension (`--by-nodes`, `--by-vdim`). It prints the time of the
setup, the mean time of `--repeat` resolutions, the time of the setup of a
second operator of the same sparsity, the number of iterations, the
relative residual and the convergence reported (the first resolution after
each setup is timed apart, since some preconditioners are built at their
first application):

```
mpirun -n 4 ./LinearSolversBenchmark -p elasticity -n 16 -c 1000 \
  -s BelosSolver -sp '{Solver=CG,RelativeTolerance=1e-8,Preconditioner={Name=MueLu}}'
```
