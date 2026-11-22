# Series updated with T1: tangent operator computed by a single call in device code (2026-11-22)

Series `78fb1d845` = series `a74d1a2ab` + T1: in device code, the tangent
operator (or prediction operator) defined by a finite strain behaviour is
computed by a single call in `computeConsistentTangentOperator`
(`computePredictionOperator`), the functions of the operators deduced from
it only converting it (`BehaviourCodeGeneratorBase.cxx`; host code
unchanged). Targeted validation and measurements of the polycrystal law
(`scratch/night-2026-11-22/validate-series3.sh`), against the series
`a74d1a2ab` (law libraries `-series2`, same law and options).

```
validation of the series 78fb1d845 Documentation of the generic-parallel interfaces
== build and install TFEL (-j8) 10:53
status=0
== targeted TFEL tests 10:53
build status=0
100% tests passed, 0 tests failed out of 6752
== sweep cuda-nvcc (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 10:55
status=0      15 OK       1 USER_CODE 
== sweep cuda-clang (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 10:56
status=0      15 OK       1 USER_CODE 
== generic-parallel cuda-nvcc -O3 10:57
behaviours           cuda-nvcc                        SUCCESS   469 s (2 passed, 0 skipped)
material-properties  cuda-nvcc                        SUCCESS    30 s (1 passed, 0 skipped)
== compile FiniteStrainSingleCrystal with openacc-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with stlpar-parunseq-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with sycl-default-acpp: status=0 1 library
== polycrystal law with the new series 11:12
law cuda float: status=0
law cuda double: status=0
law openmp float: status=0
law openmp double: status=0
== static resources of the law kernel (H100 sm_90, MI300A gfx942) 11:13
-series2 float gfx942: VGPR  128  AGPR    0  scratch   7848 B/lane  waves/SIMD 4  VGPR spill 15
-series2 float sm90: regs  255  stack   4152 B  spill st/ld 112/120 B
-series3 float gfx942: VGPR  128  AGPR    0  scratch   6040 B/lane  waves/SIMD 4  VGPR spill 19
-series3 float sm90: regs  255  stack   3912 B  spill st/ld 4/4 B
-series2 double gfx942: VGPR  128  AGPR    0  scratch  11504 B/lane  waves/SIMD 4  VGPR spill 18
-series2 double sm90: regs  255  stack   9984 B  spill st/ld 1964/2832 B
-series3 double gfx942: VGPR  128  AGPR    0  scratch   8672 B/lane  waves/SIMD 4  VGPR spill 19
-series3 double sm90: regs  255  stack   9488 B  spill st/ld 1956/2832 B
== MI300A: series with and without #1501 (overlay of LUDecomp.ixx) 11:13
-series3-no1501 float gfx942: VGPR  128  AGPR    0  scratch   6048 B/lane  waves/SIMD 4  VGPR spill 18
-series3-no1501 double gfx942: VGPR  128  AGPR    0  scratch   8728 B/lane  waves/SIMD 4  VGPR spill 17
== bench binaries 11:13
bench-gpu-double
bench-gpu-float
bench-host-double
bench-host-float
== bit-identity (last step outputs, base = -series2) 11:13
gpu float: identical
cpu float: identical
gpu double: identical
cpu double: identical
== timings, interleaved rounds (median of the medians of 'next ones') 11:13
gpu float n=1048576 -series2: 5.172 ms (5.189 5.166 5.170 5.172 5.176)
gpu float n=1048576 -series3: 5.166 ms (5.149 5.159 5.166 5.168 5.170)
cpu1 float n=8192 -series2: 31.343 ms (31.343 31.201 31.488 31.092 31.397)
cpu1 float n=8192 -series3: 31.362 ms (31.018 31.373 30.995 31.392 31.362)
gpu double n=262144 -series2: 20.499 ms (20.511 20.498 20.499 20.499 20.491)
gpu double n=262144 -series3: 20.492 ms (20.488 20.494 20.492 20.488 20.494)
cpu1 double n=8192 -series2: 29.990 ms (29.900 30.046 29.981 29.990 30.188)
cpu1 double n=8192 -series3: 30.010 ms (30.010 29.938 30.129 30.174 29.959)
== end 11:14
```

With T1, #1501 no longer grows the scratch frame of MI300A (6040 B with
it, 6048 B without in float; 8672 B and 8728 B in double): the growth
came from the AMDGPU limit on the number of basic blocks of an inlined
function (`-amdgpu-inline-max-bb`, 1100), crossed by the 12 copies of the
tangent operator of the behaviour inlined in the integration, which T1
reduces to one. The "VGPR spill" count of the remarks only covers the
kernel entry function.

Static numbers of the kernel from the analysis of T1 (graph base, base / T1):
instructions -44 % (RTX float) to -77 % (MI300A), kernel compilation time
-20 % to -70 %. RTX run time unchanged (float and double).
