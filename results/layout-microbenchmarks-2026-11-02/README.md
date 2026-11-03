# Layout of the data at the integration points: micro benchmarks (02/11/2026)

Purpose: decide how the layout of the data at the integration points (point
by point, `AoS`, or component by component, `SoA`) should be chosen, before
writing it in mfem-mgis (see the entry of 03/11/2026 of
[JOURNAL.md](../../JOURNAL.md)). The kernels are written for this purpose:
they reproduce the access patterns of the kernels of mfem-mgis, they are
not the kernels of mfem-mgis.

Machine: RTX 5080 (sm_120, nvcc 13.0, copy bandwidth 810 GB/s), CPU with 8
cores (g++ 13.3, `-O3 -march=native`, threads bound to the cores 0 to 7).
Double precision, 2 million points.

## Files

| Sources | What is measured | Results |
|---|---|---|
| `kernels.hxx`, `gpu.cu`, `cpu.cxx` | accesses only: `y = K g` at each point, no other computation | `results-gpu-5080.txt`, `results-cpu-1-core.txt`, `results-cpu-8-cores.txt` |
| `mfaction.cu`, `cpu_mf.cxx` | action of the tangent operator without assembly, hexahedra Q1: gradient of the increment from the 8 nodes, product by K, then the transposed operator B | `results-action-gpu-5080.txt`, `results-action-cpu.txt` |
| `mfaction_rt.cu`, `cpu_mf_rt.cxx` | the same, layout known at compile time or strides read at run time | `results-runtime-strides-gpu.txt`, `results-runtime-strides-cpu.txt` |

Compilation: `nvcc -O3 -std=c++20 -arch=sm_120 -DGV=9 -DNQV=8 mfaction_rt.cu`
(`GV`: 6 or 9 components of the gradients, `NQV`: points per element) and
`g++ -O3 -march=native -std=c++20 -fopenmp -DGV=9 -DNQV=8 cpu_mf_rt.cxx`,
run with `OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close taskset -c 0-7`.

## Action of the tangent operator, GPU (GB/s, strides read at run time)

| K | points per element | one thread per element, AoS (current kernels) | one thread per element, SoA | one thread per point, AoS | one thread per point, SoA |
|---|---|---|---|---|---|
| 36 components (small strain) | 8 | 593 | 569 | 591 | 599 |
| 36 components | 27 | 496 | 172 | 506 | 511 |
| 81 components (finite strain) | 8 | 339 | 255 | 430 | 723 |
| 81 components | 27 | 326 | 156 | 416 | 692 |

"One thread per point": the contribution of each point is stored in shared
memory, then summed in the order of the points by one thread per value of
the element (deterministic reduction, as the partial assembly kernels of
MFEM). All the variants give the same results to the last bit.

## Action of the tangent operator, CPU (time SoA / time AoS)

| K | points per element | 1 core | 8 cores |
|---|---|---|---|
| 36 | 8 | 1.05 (2.20 if the number of points is a power of two) | 1.29 to 1.38 |
| 36 | 27 | 1.06 to 1.10 | 1.06 to 1.07 |
| 81 | 8 | 1.95 (4.88 if the number of points is a power of two) | 1.10 (1.44) |
| 81 | 27 | 1.92 to 2.56 | 1.14 to 1.16 |

## Conclusions

1. On the GPU, the component by component layout only pays with one thread
   per point: with one thread per element, the mapping of the current
   kernels of the finite elements, it is at best equal, otherwise 1.3 to
   2.9 times slower. With one thread per point, it is twice as fast as the
   current kernel for K with 81 components, and equal for 36.
2. One thread per point is the best mapping on the GPU in both layouts
   (+27 % in AoS for K with 81 components).
3. On the CPU, the point by point layout is always better or equal. The
   component by component layout costs 5 to 45 % on 8 cores, and up to a
   factor 5 on one core when the number of points is a power of two: all
   the columns fall in the same sets of the caches.
4. Strides read at run time cost nothing measurable: no difference in SoA
   on the GPU, -23 % to +9 % in AoS without trend (scheduling of the
   compiler), ±7 % on the CPU. One version of each kernel is enough.
5. Columns not aligned on 256 bytes (2000008 points): no effect on this GPU.
6. The results do not depend on the layout or on the mapping, as long as
   the sums keep the order of the points.

The real kernels of mfem-mgis remain to be measured in both layouts.
