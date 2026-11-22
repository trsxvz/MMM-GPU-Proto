# Series updated with thelfer/tfel#1489, #1491, #1493, #1495, #1497 (night of 2026-11-22)

Targeted validation and measurements of the polycrystal law (`scratch/night-2026-11-22/validate-series.sh`), against the series before the night (law libraries built the day before, same law and options).

```
validation of the series ad6e92787 Documentation of the generic-parallel interfaces
== build and install TFEL (-j8) 04:23
status=0
== targeted TFEL tests 04:24
build status=2
98% tests passed, 12 tests failed out of 525
The following tests FAILED:
	1393 - generic-UninitializedVariablesTest_UpWard_mtest (Failed)
	1394 - generic-UninitializedVariablesTest_DownWard_mtest (Failed)
	1395 - generic-UninitializedVariablesTest_ToNearest_mtest (Failed)
	1396 - generic-UninitializedVariablesTest_TowardZero_mtest (Failed)
	7624 - tdls-TDLSPartialJacobianInvertTest_UpWard_mtest (Failed)
	7625 - tdls-TDLSPartialJacobianInvertTest_DownWard_mtest (Failed)
	7626 - tdls-TDLSPartialJacobianInvertTest_ToNearest_mtest (Failed)
	7627 - tdls-TDLSPartialJacobianInvertTest_TowardZero_mtest (Failed)
	7628 - tdls-TDLSPartialJacobianInvertTest2_UpWard_mtest (Failed)
	7629 - tdls-TDLSPartialJacobianInvertTest2_DownWard_mtest (Failed)
	7630 - tdls-TDLSPartialJacobianInvertTest2_ToNearest_mtest (Failed)
	7631 - tdls-TDLSPartialJacobianInvertTest2_TowardZero_mtest (Failed)
Errors while running CTest
Output from these tests are in: /home/t/Documents/GitHub/mfem-mgis-gpu-proto/build/tfel-mmm-gpu/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
== sweep cuda-nvcc (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 04:24
status=0      15 OK       1 USER_CODE 
== sweep cuda-clang (FiniteStrainSingleCrystal|SaintVenantKirchhoff) 04:25
status=0      15 OK       1 USER_CODE 
== generic-parallel cuda-nvcc -O3 04:26
behaviours           cuda-nvcc                        SUCCESS   492 s (2 passed, 0 skipped)
material-properties  cuda-nvcc                        SUCCESS    32 s (1 passed, 0 skipped)
== compile FiniteStrainSingleCrystal with openacc-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with stlpar-parunseq-nvhpc-gpu: status=0 1 library
== compile FiniteStrainSingleCrystal with sycl-default-acpp: status=0 1 library
== polycrystal law with the new series 04:41
law cuda float: status=0
law cuda double: status=0
law openmp float: status=0
law openmp double: status=0
== bench binaries 04:42
bench-gpu-double
bench-gpu-float
bench-host-double
bench-host-float
== bit-identity (last step outputs, base = series before tonight) 04:42
gpu float: identical
cpu float: different
gpu double: identical
cpu double: different
== timings, interleaved rounds (median of the medians of 'next ones') 04:42
gpu float n=1048576 -e1-base: 7.434 ms (7.452 7.414 7.406 7.434 7.441)
gpu float n=1048576 -night-series: 6.895 ms (6.884 6.895 6.889 6.905 6.905)
cpu1 float n=8192 -e1-base: 44.358 ms (45.096 44.543 43.880 44.358 43.790)
cpu1 float n=8192 -night-series: 42.263 ms (41.421 44.110 42.263 41.906 42.581)
gpu double n=262144 -e1-base: 21.983 ms (21.983 21.990 21.975 21.979 22.001)
gpu double n=262144 -night-series: 21.112 ms (21.125 21.100 21.106 21.126 21.112)
cpu1 double n=8192 -e1-base: 41.697 ms (42.248 41.650 41.777 40.266 41.697)
cpu1 double n=8192 -night-series: 38.412 ms (38.634 38.411 37.043 38.412 38.710)
== end 04:42
```

Targeted tests rerun after building the right targets (the first run stopped on a wrong target name, the test libraries were stale): 100 % of 6748 tests passed (generic-interface behaviours in the four rounding modes, TDLS behaviours, conversions, sparse products, TinyMatrixSolve, fixed size arrays, Newton-Raphson, mfront-DSLTest).

CPU outputs (law compiled with -march=native, hence FMA): float, 76 of 737280 values differ, at most 1.8e-12 of the largest value; double, 1.9 % of the values differ, at most 5.5e-17 of the largest value (contractions of the sparse products). GPU outputs identical to the bit, float and double.
