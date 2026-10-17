# Solvers of PETSc with MFEM on a GPU: wrong solutions

The example `examples/petsc/ex1p.cpp` of MFEM 4.10 (unmodified copy, as
the options file `rc_ex1p_device` and the mesh `star.mesh`), whose header
gives the device runs `mpirun -np 4 ex1p -pa -d cuda --petscopts
rc_ex1p_device`, solves the Poisson problem with the conjugate gradient of
PETSc. On the GPU, the conjugate gradient stops after a few iterations
with `CONVERGED_ATOL`, its residual being seen as zero, and the solution
is wrong (zero on one process). On the CPU, the same problem is solved
normally. The files `ex1p.cpp` and `linalg/petsc.cpp` are the same in the
main branch of MFEM (10/2026).

## Environment

- MFEM 4.10.0 built with CUDA 13.0.2 (`-arch=sm_120`), MPI, hypre 3.2.0
  and PETSc (spack, `share/mfem/config.mk` of the installation);
- PETSc 3.26.0 built with CUDA and Kokkos (`PETSC_HAVE_CUDA`,
  `PETSC_HAVE_DEVICE`, double precision);
- OpenMPI 5.0.10, one RTX 5080 (16 GB), driver CUDA 13.2.

## Run

```
reproducers/petsc-gpu-vectors/run.sh [MFEM installation]
```

The script builds `ex1p` with the compiler and the flags of the
installation of MFEM (`makefile`, which includes its `config.mk`), in
`scratch/reproducers/petsc-gpu-vectors`, then runs it on 1 and 2
processes with `-d cpu` and `-d cuda` (and `-pa -d cuda`: operator applied
without assembly, given to PETSc as a shell matrix), with the options of
`rc_ex1p_device` (conjugate gradient, default preconditioner: block Jacobi
with ILU) and `-ksp_monitor_true_residual`. It prints the first residual,
the reason of the end of the resolution and the largest value of the
solution written by the example (`sol.000000`, first process).

## Results (17/10/2026)

| Run | First residual | End of the resolution | Largest value of the solution |
|---|---|---|---|
| 1 process, cpu | 2.13e-2 | `DIVERGED_ITS` after 200 iterations | 0.3226 |
| 1 process, cuda | 0 | `CONVERGED_ATOL` after 0 iteration | 0 |
| 2 processes, cpu | 2.13e-2 | `DIVERGED_ITS` after 200 iterations | 0.3216 |
| 2 processes, cuda | 2.13e-2 | `CONVERGED_ATOL` after 87 iterations | 1.27e-4 |
| 2 processes, cuda, `-pa` | 9.88e-3 | `CONVERGED_ATOL` after 3 iterations | 5.80e-5 |

On one process with the GPU, `-ksp_monitor_true_residual` gives at the
iteration 0 a preconditioned (natural) residual of 0 and a true residual
of 1.65e-2, the same as on the CPU: the right hand side is seen, the
preconditioned residual is zero. (The 200 iterations on the CPU are the
limit of `rc_ex1p_device`, which does not choose the preconditioner: the
residual decreases normally, down to 2.2e-6 relatively.) The program
ends normally in all the runs.

## Leads

- MFEM wraps its vectors in the vectors of PETSc without copy
  (`PetscParVector`, `linalg/petsc.cpp`) and sets the internal offload mask
  of the vectors of PETSc itself, through the private header
  `petsc/private/vecimpl.h` (`__mfem_VecSetOffloadMask`, l. 5512-5519 of
  `linalg/petsc.cpp`, and its callers). A change of the handling of the
  device arrays of the vectors of PETSc (`VecCUDA`, `VecCUPM`) would make
  PETSc read an array which is not up to date (the array of the host while
  the values are on the device, or the converse).
- To check: the same problem with PETSc alone (`-vec_type cuda -mat_type
  aijcusparse`, for instance a tutorial of KSP of PETSc) to rule out PETSc;
  the offload mask and the arrays of the host and of the device of the
  vectors given to the preconditioner (`PCApply`) at the iteration 0; the
  public functions of PETSc (`VecCUDAPlaceArray`, `VecCUDAResetArray`,
  `VecCUDAGetArray`...) instead of the private field; an older version of
  PETSc, to know whether the behaviour changed.
- In mfem-mgis, the solvers of PETSc compute on the host since 15/10/2026
  (`src/PETScLinearSolvers.cxx`): the matrices are copied in the AIJ format
  of the host and the vectors through the arrays of the host. The same
  vectors seen as zero were met there, with also an abort of the program in
  `PetscFinalize` when the memory pool of PETSc on the GPU was freed: this
  abort does not happen with `ex1p`.
