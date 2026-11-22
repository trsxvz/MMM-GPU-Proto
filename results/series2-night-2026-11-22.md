# Series updated with the runtime, LU and tangent conversion changes (morning of 2026-11-22)

Series `a74d1a2ab` (tip `709393517` when measured, same tree except the
release notes) = series `ad6e92787` + runtime commits 1-3 of the
generic-parallel interface + thelfer/tfel#1500 and #1501 (LU) + conversions
of the finite strain tangent operators (generator, 62 edges, breadth first
search of the path, wrappers of the generic interface). Targeted validation
and measurements of the polycrystal law
(`scratch/night-2026-11-22/validate-series2.sh`), against the series
`ad6e92787` (law libraries `-night-series`, same law and options).

```
validation of the series 709393517 Documentation of the generic-parallel interfaces
== build and install TFEL (-j8) 09:10
status=0
== targeted TFEL tests 09:11
build status=0
99% tests passed, 1 tests failed out of 6752
The following tests FAILED:
	7756 - FiniteStrainBehaviourTangentOperator17 (Not Run)
Errors while running CTest
Output from these tests are in: /home/t/Documents/GitHub/mfem-mgis-gpu-proto/build/tfel-mmm-gpu/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
== sweep cuda-nvcc (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 09:13
status=0      15 OK       1 USER_CODE 
== sweep cuda-clang (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 09:15
status=0      15 OK       1 USER_CODE 
== generic-parallel cuda-nvcc -O3 09:15
behaviours           cuda-nvcc                        SUCCESS   509 s (2 passed, 0 skipped)
material-properties  cuda-nvcc                        SUCCESS    30 s (1 passed, 0 skipped)
== compile FiniteStrainSingleCrystal with openacc-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with stlpar-parunseq-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with sycl-default-acpp: status=0 1 library
== polycrystal law with the new series 09:32
law cuda float: status=0
law cuda double: status=0
law openmp float: status=0
law openmp double: status=0
== static resources of the law kernel (H100 sm_90, MI300A gfx942) 09:32
-night-series float gfx942: VGPR  128  AGPR    0  scratch   6704 B/lane  waves/SIMD 4  VGPR spill 10
-night-series float sm90: regs  255  stack   4432 B  spill st/ld 272/500 B
-series2 float gfx942: VGPR  128  AGPR    0  scratch   7848 B/lane  waves/SIMD 4  VGPR spill 15
-series2 float sm90: regs  255  stack   4152 B  spill st/ld 112/120 B
-night-series double gfx942: VGPR  128  AGPR    0  scratch   9128 B/lane  waves/SIMD 4  VGPR spill 13
-night-series double sm90: regs  255  stack  10328 B  spill st/ld 3108/3404 B
-series2 double gfx942: VGPR  128  AGPR    0  scratch  11504 B/lane  waves/SIMD 4  VGPR spill 18
-series2 double sm90: regs  255  stack   9984 B  spill st/ld 1964/2832 B
== bench binaries 09:33
bench-gpu-double
bench-gpu-float
bench-host-double
bench-host-float
== bit-identity (last step outputs, base = series ad6e92787) 09:33
gpu float: identical
cpu float: identical
gpu double: identical
cpu double: identical
== timings, interleaved rounds (median of the medians of 'next ones') 09:37
gpu float n=1048576 -night-series: 6.866 ms (6.864 6.866 6.866 6.877 6.862)
gpu float n=1048576 -series2: 5.176 ms (5.186 5.169 5.174 5.177 5.176)
cpu1 float n=8192 -night-series: 31.732 ms (31.732 31.642 32.087 31.644 32.499)
cpu1 float n=8192 -series2: 31.101 ms (31.101 31.015 31.106 31.375 30.995)
gpu double n=262144 -night-series: 21.095 ms (21.078 21.101 21.095 21.095 21.104)
gpu double n=262144 -series2: 20.503 ms (20.506 20.503 20.499 20.500 20.504)
cpu1 double n=8192 -night-series: 30.946 ms (31.044 30.993 30.946 30.847 30.868)
cpu1 double n=8192 -series2: 30.933 ms (30.933 31.029 30.917 31.175 30.822)
== end 09:42
```

`FiniteStrainBehaviourTangentOperator17` was missing from the list of the
targets built by the script (fixed since): built and run afterwards, the 17
tests of the conversions of the tangent operators pass.

Contribution of #1501 alone (series with the `LUDecomp.ixx` of `ad6e92787`
given by a header overlay, `-series2-no1501`):

| | series | without #1501 |
|---|---|---|
| RTX 5080 GPU float (3 rounds) | 5.185-5.196 ms | 5.782-5.823 ms |
| RTX 5080 GPU double (3 rounds) | 20.509-20.518 ms | 20.654-20.662 ms |
| H100 float: stack, spills st/ld | 4152 B, 112/120 B | 4144 B, 108/116 B |
| H100 double: stack, spills st/ld | 9984 B, 1964/2832 B | 9992 B, 2356/3112 B |
| MI300A float: scratch, VGPR spills | 7848 B, 15 | 6424 B, 15 |
| MI300A double: scratch, VGPR spills | 11504 B, 18 | 8936 B, 18 |

The scratch frame of MI300A grows with #1501 (a large calling function
is no longer inlined by clang; the static scratch traffic decreases by
14 to 24 % according to the analysis of the generated code); to be
measured on MI300A, analysis in progress.
