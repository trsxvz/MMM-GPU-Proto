# Pitfalls

This page lists the pitfalls met during the work, each with its symptom
and the workaround adopted. Compiler defects are worked around in a
targeted way, never by removing a GPU feature, and each one has a
reproducer in `results/nvhpc-bugs/` or `reproducers/`. The sources
are `JOURNAL.md`, `results/nvhpc-bugs/README.md`,
`reproducers/README.md`, `deps/spack/README.md`, `docs/gpu-solvers.md`,
`docs/plan-precision.md` and the header of `scripts/build-hip-amd.sh`.

A general lesson: several of these defects give **wrong results without
any message**. Only numerical comparisons against the `generic`
interface on the CPU, register counts of the kernels, and tests checked
against deliberately wrong code revealed them.

## nvcc (CUDA 13)

- **Calls to host functions from device code.** Symptom: nvcc only warns
  (#20011, #20014) and removes the corresponding code path: empty kernel
  or wrong result. Workaround: `-diag-error 20011,20014` in the
  configurations that use nvcc.
- **`std::string_view` built from a C string in device code** (with
  libstdc++ 14). Symptom: the whole kernel is emptied, without error or
  warning (Norton kernels had 4 registers and computed nothing).
  Workaround: `const char*` in the bound checks; a small string type
  without `std::char_traits` works with all six compilers
  (`reproducers/device-strings`).
- **`std::fpclassify` in device code.** Symptom: empty kernel.
  Workaround: `tfel::math::ieee754::fpclassify`.
- **Functions using `__builtin_constant_p`.** Symptom: compiled without
  diagnostic into empty functions in device code, so `std::advance`,
  `std::next` and `std::prev` do nothing on the GPU (cubic spline
  interpolation gives NaN). Workaround: internal helpers in
  `CubicSpline.ixx`.
- **Local `static constexpr` objects.** Symptom: the device copy of such
  an object whose constructor is called is filled with zeros, without
  diagnostic. Workaround: namespace scope `__device__ constexpr`
  variables (`tfel::utilities::getConstantObject`).
- **`auto i = decltype(s){};`** Symptom: cudafe++ rewrites it into an
  invalid expression for the host compiler. Workaround: `size_type{}`.
- **NVVM optimiser.** Symptom: the stress written by an initialize
  function of a finite strain behaviour was lost when the end of step
  state was copied after the construction of the point view; correct
  with `-G`, `-Xcicc -O0` and clang CUDA, and changing with added
  `printf`. Workaround: the copy is made in the constructor of
  `PointView`.
- **Kernel launcher.** Symptom: nvcc cannot generate the host launcher of
  a kernel whose non-type template parameter is a class. Workaround:
  sizes passed through a type (`PointTraits`).
- **Anonymous C structures.** Symptom: cudafe crashed on the destructor
  of a `std::vector` of anonymous structures. Workaround: named
  structures in `BatchData.h`.
- **Extended lambdas.** Symptom: nvcc forbids initialisers in the
  captures of `__host__ __device__` lambdas. Workaround: variables
  declared before the kernels and captured by name.
- The defects of the local `static constexpr` objects and of
  `decltype(s){}` were still present from nvcc 13.0.88 to 13.4.92.

## nvc++ (NVIDIA HPC SDK 25.11)

- **Signed loop counters.** Symptom: on the GPU at `-O2` and `-O3`, a
  loop `for (int k = 0; k != N; ++k)` runs N+1 times. Workaround:
  unsigned counters in TFEL and in the generated code (user code of
  `.mfront` files stays exposed).
- **Unified memory mode** (default on a machine with HMM). Symptom: a
  global store to a stack address in the constructor of the behaviour.
  Workaround: `-gpu=mem:managed`.
- **Default initialiser `= {}` of an array member.** Symptom: invalid LLVM
  constant. Workaround: member initialised in the constructor (fixed in
  26.9).
- **Calls through a pointer to member function.** Symptom: silently
  ignored on the GPU, at all optimisation levels. Workaround: MFront
  generates named function objects.
- **`std::variant` visit.** Symptom: "Unexpected branch type", invalid
  LLVM aliases. Workaround: `TrivialGenTypeBase` for the tangent
  operator.
- **Link of the device code (`-stdpar=gpu`).** Symptom: "nvlink error:
  Multiple definition of `__staticinit_...`", appearing at random in a
  library of several behaviours. Workaround: `-gpu=nordc`.
- **OpenMP offloading (`-mp=gpu`).** Several defects, each with a
  workaround:
  - the consistent tangent of implicit DSLs returned false: a local
    object held a `bool&` to a local flag, each object now owns its flag;
  - a static block holding the text of an `assert` was not found at load
    time: the redundant `assert` is removed;
  - the zero initialisation of an array member removed the kernel code:
    a loop in the constructor replaces it;
  - a lambda capturing a `std::array` of pointers read wrong pointers:
    the pointers are captured one by one;
  - a nested lambda capturing by reference did not write: it captures by
    copy;
  - an internal error "Computed goto control expression must be a 32-bit
    integer": the tensor view is written without `switch`.
- **Global and static variables.** `-stdpar=gpu` forbids them in GPU
  code (explicit error): `getConstantObject` returns a copy.
- **Environment.** `-mp=gpu -gpu=mem:managed` does not link
  `libnvhpcmanaux`, so a program that frees managed memory from g++ code
  can stop on `free(): invalid pointer` (link with `-lnvhpcmanaux`).
  `libnvcpumath.so` uses `sincos` from `libm` without depending on it.
- **Release notes.** Those of NVHPC 26.1 to 26.9 list no fixed defects:
  only the reproducers tell. Of the first five reproducers, only the
  array initialiser one is fixed in 26.9.
- **In mfem-mgis**, the `stdpar` and OpenMP offload libraries of nvc++
  fail: the OpenMP runtime of NVHPC conflicts with the `libgomp` of MFEM
  and mfem-mgis compiled by GCC.

## clang CUDA, HIP and SYCL

- **Attribute order.** Symptom: clang CUDA compiled no behaviour, since
  `TFEL_HOST_DEVICE [[nodiscard]]` is refused (a C++11 attribute must
  come first). Workaround: 285 corrections in 41 headers and in the
  generators.
- **Inconsistent annotations.** Symptom: clang refuses a definition
  without `__host__ __device__` whose declaration has them; nvcc accepts
  it. Workaround: 85 definitions corrected, and every header of
  `TFEL/Math` and `TFEL/Material` compiled alone with clang CUDA
  (`scripts/clang-cuda-header-scan.sh`).
- **CUDA versions.** clang 20 only knows CUDA up to 12.8 and fails with
  CUDA 13.0. Workaround: the CUDA 12.9 of the HPC SDK, chosen by the
  CMake variable `TFEL_CLANG_CUDA_PATH`, with its libraries after those
  of CUDA 13 in `LD_LIBRARY_PATH`.
- **hipcc on the NVIDIA platform** calls nvcc and defines neither
  `__HIP__` nor `__HIPCC__`, so the code silently took the CUDA path.
  Workaround: include `hip/hip_runtime.h` first. hipcc also does not give
  the HIP headers to nvcc: older system headers (HIP 5.7) were taken
  instead of those of HIP 7.2. Workaround: the HIP include directory is
  added to the hipcc options when it is not an implicit directory of the
  C++ compiler.
- **MFront and `HIP_COMPILER`.** MFront reads this variable to choose the
  HIP compiler; `HIP_COMPILER=nvcc`, meant for hipcc, made it call nvcc.
  Workaround: variable removed from the environment.
- **AdaptiveCpp.** Its generic compiler is single pass, like nvc++:
  `TFEL_DEVICE_CODE` is defined in the whole translation unit. Virtual
  functions require RTTI, absent from device code ("Unresolved extern
  variable" from ptxas): `MechanicalBehaviour` lost its virtual
  functions. Reductions with barriers fail on large kernels run on the
  CPU: SYCL uses atomic updates instead. A recursion in TFEL made one of
  its passes loop forever.
- **Kokkos `nvcc_wrapper`.** `-I dir` with a space is taken for a source
  file (use `-Idir`); a doubled `-std` prints a warning on standard
  output that corrupts dependency files.

## Host compilers and floating point

- **`std::for_each(par_unseq, ...)` over `std::views::iota`** runs on a
  single thread with libstdc++ (the iterator is not random access in the
  C++17 sense). Workaround: a counting iterator.
- **g++ at `-O2`** vectorises only with a very restrictive cost model:
  the polycrystal behaviour was 17 % slower than with the `generic`
  interface. Workaround: the flags of `tfel-config --oflags`.
- **FMA contraction.** Different compilers contract differently, and the
  local Newton solvers amplify it. A 66 % gap of a Hill behaviour on the
  GPU also appeared on the CPU with `-ffp-contract=fast`: an ill
  conditioned problem, not a GPU defect. Rule: the `generic` and
  `generic-parallel` libraries compiled with FMA on the CPU must agree
  before any GPU comparison (`scripts/roundoff-check.sh`,
  `scripts/roundoff-triage.py`). In float on the GPU, see the finding of
  14/10/2026 in [Corrections](corrections.md).

## GPU runtime and memory

- **GPU thread stack.** CUDA threads have 1 KiB of stack by default, not
  increased automatically, and the static estimate ignores internal
  routines such as double division (stack overflow seen with
  compute-sanitizer). Workaround: the backends compute the size from the
  kernel attributes and raise the limit (`GPUStackSize.hxx`,
  `DeviceStackSize.hxx`).
- **Large stacks.** A polycrystal behaviour needs 113 to 121 KiB per
  thread, so 14.8 to 15.9 GB are reserved on the 16 GB GPU whatever the
  number of points: another process then causes an "out of memory". Such
  a behaviour must have the GPU to itself.
- **Statically linked CUDA runtime.** Each library has its own instance;
  a shared function taking the kernel as an argument was resolved to one
  library, and the kernels of the others were unknown ("invalid resource
  handle"). Workaround: the kernel is a template parameter.
- **TBB threads at exit.** Libraries linked to TBB were unloaded while its
  threads were still running (intermittent crash at the end of the
  process). Workaround: `-Wl,-z,nodelete`.
- **Kokkos** must be a shared library, compiled without relocatable
  device code, shared by the caller and the generated libraries. In
  mfem-mgis, linking the static library to Kokkos leaked the Kokkos link
  options to all client programs; the Kokkos code is now in a small
  separate shared library.
- **Heterogeneous memory management (HMM).** The NVIDIA driver of the
  development machine lets the GPU read ordinary host memory, which hides
  location errors (an array of argument pointers read on the caller's
  stack, a test writing host memory from a kernel). Workaround: final
  validations with HMM disabled (`uvm_disable_hmm=1`).
- **Validity of MFEM memory.** MGIS writes through raw pointers that MFEM
  does not see. Every host access goes through a synchronisation first,
  and pointers given to MGIS are obtained by `Read`, `Write` or
  `ReadWrite` just before each call. The `debug` device of MFEM, which
  runs the kernels on the host and protects the host copies, detects
  missing synchronisations.
- **Stale host copies.** The `operator[]` of an MFEM vector reads the
  host copy without making it valid (L2 projections gave zeros on the
  GPU); `MakeTRef` creates an alias with frozen validity flags (the
  ParaView export overwrote the unknowns); Dirichlet conditions wrote the
  unknowns on the host without declaring it. Workaround: `HostRead` and
  `HostReadWrite`.
- **Umpire for host memory.** With MFEM compiled with Umpire, host memory
  is also allocated by Umpire once a device is configured, and every
  temporary of the host element loops goes through its allocation table:
  the bubble case took 68 s on the GPU against 31 s on one CPU core.
  Workaround: standard host memory, Umpire kept for the GPU.
- **hypre memory pool.** hypre 3.2 reserves 4 GiB on the GPU at once, per
  process. Workaround: option `--hypre-device-pool-size` of mfem-mgis.

## MFEM, hypre and solver libraries

- **`mfem::forall`** only launches a CUDA kernel in a file compiled by
  nvcc; otherwise the body runs on the host. Files with kernels are
  compiled as CUDA (or HIP).
- **Partial assembly of MFEM 4.10** silently ignores boundary
  integrators and stops on domain integrators with markers; elements are
  in lexicographic order. mfem-mgis uses its own operator.
- **`QuadratureInterpolator::MultTranspose`** is not implemented in MFEM
  4.10: the internal forces kernel is written in mfem-mgis.
- **Geometric factors.** `Mesh::GetGeometricFactors` gives nodes to the
  mesh, which made the remaining element loops 50 % slower; the shape
  function cache of the global transformation element is indexed by the
  address of the integration rule, so two rules allocated at the same
  address gave wrong Jacobians. Workaround: Jacobians computed by a
  kernel from the vertices or nodes.
- **`HypreParMatrix` block diagonal constructor** moves the diagonal entry
  to the head of each row in place, without touching the sorted flag:
  a pattern with sorted columns was wrong in parallel. The pattern now
  follows the order of hypre.
- **`EliminateRowsCols`** copies the matrix to the host (122 ms per call
  on the beam). Replaced by `EliminateBC` (`DIAG_ONE`), in place on the
  GPU.
- **BoomerAMG on the GPU.** Its defaults depend on the device at
  construction, so the AMG must be created after `mfem::Device`; the
  `Elasticity` strategy stops the program (`System` used instead).
- **`HypreILU`** is recreated with the default fill level at each new
  operator; with hypre 3.2.0, ILU(k) with k at least 1 crashes on the GPU
  (`reproducers/hypre-ilu-gpu`). Workaround: settings reapplied, level 0
  on the GPU; ILUT works.
- **ParaSails** only computes on the host and crashes when hypre computes
  on the GPU: refused there.
- **hypre in single precision** calls the double precision LAPACK
  routines `dpotrf_` and `dpotrs_`: ParaSails gives NaN. The ParaSails
  test stays in double only.
- **cuDSS.** MFEM chooses the reuse of the analysis only before the first
  matrix, stops on a different size, and only checks the number of
  coefficients; MUMPS reuse does not check the pattern at all.
  Workaround: mfem-mgis compares the pattern and recreates the solver.
  cuDSS does not report factorisation failures, so each solution is
  checked by its residual. The distributed analysis of cuDSS 0.7.1 fails
  for some matrices that MUMPS solves (`reproducers/cudss-distributed`).
- **AmgX.** Resources hold global objects freed by the destruction of any
  resources: one resource per process. `exception_handling=1` makes AmgX
  abort. `AMGX_vector_set_zero` resizes the vector without the halo of a
  distributed matrix. Convergence has floors (1e-20 absolute, 1e-12
  relative in double, 1e-6 in single), handled by scaling the right hand
  side by a power of 2. `libmfem.so` contains its own unpatched static
  AmgX: the patched `libamgxsh.so` must be linked before MFEM, and
  solvers destroyed before `MPI_Finalize`. The mixed precision mode
  `dDFI` is no longer supported by AmgX 2.5.

## MPI and process placement

- **Thread binding.** `scripts/env.sh` exports `OMP_PROC_BIND=spread`
  and `OMP_PLACES=threads` (for Kokkos and AdaptiveCpp). libgomp then
  pins the main thread to the first processor: the `stlpar` backend saw a
  single processor, and `openmp` with 8 threads used 4 cores. With these
  variables, `mpirun --bind-to none taskset -c ...` also puts all
  processes on the same core. Workaround: measurement scripts remove
  them.
- **Several processes on one GPU.** Resolution 9 to 14 times slower with
  two processes on the same GPU (not analysed). One process per GPU.
- **Errors on a single process.** An exception on one process left the
  others waiting in a collective operation. Workaround: failures made
  collective before the next collective operations
  (`CollectiveFailure`).
- **OpenMPI with CUDA.** Its internal hwloc enabled the OpenCL plugin
  with an invalid path. Workaround: `--disable-opencl` in the recipe.
  hypre is built without `+gpu-aware-mpi`: its exchanges go through the
  host.
- **OpenMPI 5.0.10 in the continuous integration.** Concurrent `mpirun`
  collide on job identifiers. Workaround: `OMPI_MCA_btl=self,tcp`.

## Build systems, spack and environment

- **Environment leaks.** The `LD_LIBRARY_PATH` of the shell contained an
  older spack view: the `mfront` of the build tree loaded an older
  `libTFELMFront`, and the first test pass tested the old MFront. Spack
  views also provided older TFEL headers, and `CMAKE_PREFIX_PATH` found
  the spack TFEL first. Workaround: `scripts/env.sh` neutralises these
  environments; builds use a clean environment.
- **Cached CMake substitutions** of TFEL: changing a value requires
  `cmake -U "tfel-substitutions*" .`.
- **Test binaries.** mfem-mgis test executables are outside `all` (target
  `build-tests`) and the library is static: without rebuilding
  `build-tests`, `ctest` reruns old binaries. The `BehaviourTest` library
  must be built explicitly. In MGIS, binding tests are built by `check`:
  a stale executable hid a mutation of the Fortran code.
- **`tfel-check`** reports as successful a test whose requirements are
  not met; the test script now counts skipped tests separately.
- **Spack.** A cache shared with other spack installations breaks the
  patch index (`SPACK_DISABLE_LOCAL_CONFIG` and an own cache); the BLAS
  of the HPC SDK was chosen by default (OpenBLAS imposed); Kokkos 5 has
  no `cuda_lambda` variant required by the PETSc recipe (recipe
  adapted); PETSc needs `+openmp` with a Kokkos built with OpenMP; MFEM
  needs `cxxstd=20` for Umpire 2026; kokkos-kernels needs Kokkos `+pic`;
  hypre refuses `+int64` with CUDA ("Use +mixedint"), so mfem-mgis is
  built with `~int64`; AmgX 2.5.0 and the AdaptiveCpp package of spack
  were too old for CUDA 13 (newer versions in the prototype recipes).
- **CUDA prefix.** With `bin` as a link, nvcc took the headers of the SDK,
  without those of the math libraries (`cublas_v2.h` not found).
  Workaround: `bin` as a directory of links
  (`deps/spack/make-nvhpc-cuda-prefix.sh`). CUDA 13.0 needs a driver 580
  or newer.
- **Compiling MFEM 4.10 with CUDA 13.** `ptxas` takes about 22 GB for
  `fem/lor/lor_batched.cpp`, and `cicc` did not finish
  `fem/integ/bilininteg_convection_ea.cpp` (more than 2 h 30 and 17 GB).
  Workarounds: a memory governor that suspends all build processes but
  the largest (`scripts/memory-governor.sh`), and a patch that limits the
  3D specialisations of that kernel to order 3. Also, `-lnvToolsExt`,
  absent from CUDA 13, was removed, and the static gslib library had to
  be relocatable.
- **AMD GPU (ROCm from its Python packages).** An Ubuntu package installs
  older HIP headers in `/usr/include/hip`, read before those of ROCm:
  ROCm headers are given explicitly. `scripts/env.sh` exports
  `HIP_PLATFORM=nvidia`, which hipcc follows: MFEM is compiled with
  `HIP_PLATFORM=amd`. Access to the integrated AMD GPU is reserved to the
  `render` group.
- **Integrated AMD GPU (gfx1036, ROCm 10.1).** The writes of a finished
  kernel to device memory are lost when the GPU stays idle about 300 µs
  before the next operation, for instance while the code of a shared
  library is loaded at the first launch of one of its kernels. In
  mfem-mgis, imposed components were not zeroed and hypre faulted.
  Workaround for the tests: kernels and copies serialised
  (`AMD_SERIALIZE_KERNEL=3`, `AMD_SERIALIZE_COPY=3`). Reproducer
  `reproducers/hip-idle-write-loss`, reported as
  [ROCm/rocm-systems#13469](https://github.com/ROCm/rocm-systems/issues/13469).

## Single precision

- **Coefficients out of range.** A behaviour whose coefficients leave the
  range of `float` cannot run in float (Norton test law with `A = 8e-67`
  in Pa gave `A = 0`). Workaround: rewritten as `A (seq/K)^n`; user
  behaviours will often have to be written in MPa.
- **Double literals.** An expression with a double literal and a float
  object is double, which TFEL refuses to assign to a float object; mixed
  `std::max` and `std::min` calls fail to deduce. Of 711 MFront files
  that compile in double, 131 failed only in float.
- **Criteria.** Tolerances written for double (1e-14) make no sense in
  float. Too loose a criterion (100 ε) declared a Norton increment
  converged before the first iteration, with a zero stress.
- **Residual floor.** In finite strain, the deformation gradient, close
  to identity, is only known within one ulp (6e-8 below 1, 1.2e-7 above),
  so stresses are only known within about 1.5e4 to 3e4 Pa for UO2: a
  residual relative to the initial residual of a step cannot be reached.
  The relative floor also grows with the problem size.
  Workaround: an absolute tolerance about ten times the floor.
- **Small strains in finite strain.** In the polycrystal of mm-opera-hpc,
  in float, the axial stress exceeds the double one by 6.1e-3 at an
  axial strain of 5e-5, and by 1.0e-3 at 2e-4; with a step of 1 s, the
  behaviour fails and the computation stops
  (`results/precision-polycrystal-2026-10-14`). The deformation gradient,
  stored in float close to identity, only moves by one ulp (1e-3 of a
  strain of 5e-5, a stress quantum of about 1.5e4 Pa). No workaround yet:
  an interface based on the displacement gradient H = F - I is a pending
  decision.
- **Fixed points below the resolution of float.** The fixed point of
  mm-opera-hpc on the lateral macroscopic stresses corrects the lateral
  components of F, just below 1, by about (1 + ν) (1 - 2 ν) S / E. In
  float, this correction is rounded to zero below a residual of about
  1.5e4 Pa, above the tolerance of 1e4 Pa: the first step of the cermet
  stagnated or oscillated, was cut in two, and the stress at 0.2 s was
  2.7 % away from double, an effect of the time step. Workaround: the
  fixed point takes at least twice its resolution, computed from F and
  the elastic moduli (30 804 Pa here, negligible in double)
  (`results/cermet-float-2026-10-14`).
