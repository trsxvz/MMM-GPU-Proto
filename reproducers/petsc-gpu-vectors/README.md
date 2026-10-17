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

## Cause

MFEM is built with Umpire: its default memory type on the device is then
`DEVICE_UMPIRE` (`host-umpire,cuda-umpire`). `PetscParVector::SetVecType_`
(l. 399 of `linalg/petsc.cpp`) and `PetscParMatrix::MakeWrapper` (l. 1382)
only test the memory types `DEVICE` and `MANAGED`: with `DEVICE_UMPIRE`,
the vectors given to PETSc are vectors of the host (`VECSTANDARD`), while
the matrix is a matrix of cuSPARSE (`MATAIJCUSPARSE`). The preconditioner
writes its result on the GPU, then the conjugate gradient reads the array
of the host, which stayed at zero. At the iteration 0, the vectors `B`,
`X`, `R` and `Z` of the conjugate gradient are of type `seq`, and z'r is
zero on the arrays of the host and 4.6e-4 on the arrays of the device.
With `MFEM_MEMORY=cuda` (memory of the device without Umpire:
`MFEM_MEMORY=cuda reproducers/petsc-gpu-vectors/run.sh`), the runs on the
GPU give the results of the runs on the CPU. The offload mask written
by MFEM, first suspected, is not the cause.

The fix of MFEM (all the memory types of the device, except the debug
one, give vectors of PETSc on the device) is proposed in
[mfem/mfem#5646](https://github.com/mfem/mfem/pull/5646).

## In mfem-mgis

Since `patches/mfem-mgis/0272-*`, the solvers and the preconditioners of
PETSc compute on the GPU when the linear solvers are placed there
(`--execution gpu`), without the vectors of PETSc of MFEM: the operator is
converted to the AIJ format of cuSPARSE and the vectors of PETSc of the GPU
wrap the memory of MFEM on the device (`VecCreateSeqCUDAWithArray`,
`VecCreateMPICUDAWithArray`, and their HIP versions). This
works with the MFEM 4.10 installed (with or without the fix above). The
abort in `PetscFinalize` met before was another defect: PETSc 3.26 keeps
its pools of the GPU (events, memory) in static variables of functions,
built after the static finalizer of mfem-mgis and thus destroyed before
it at the exit; PETSc is now finalized by a function registered with
`std::atexit` after the first resolution on the GPU.
