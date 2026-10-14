# Writes of a kernel lost after a short idle time of the GPU

ROCm 10.1 (HIP 7.16, Python packages of TheRock), integrated GPU gfx1036
(Ryzen 7 9800X3D), Ubuntu 24.04 with Linux 7.0. Reported to AMD:
[ROCm/rocm-systems#13469](https://github.com/ROCm/rocm-systems/issues/13469).

The writes of a finished kernel to device memory are lost when the GPU
stays idle for about 300 µs before the next operation. A later copy reads
the old content, even after `hipDeviceSynchronize`. It looks like the
modified lines of the L2 cache are dropped, without a write back, when the
GPU powers down on idle.

    make
    sg render -c "./idle-loss 100"        # 0 / 20 trials lost
    sg render -c "./idle-loss 1000"       # 20 / 20 trials lost
    sg render -c "./idle-loss 1000 sync"  # synchronised after the launch: 0 / 20
    sg render -c "AMD_SERIALIZE_KERNEL=3 ./idle-loss 1000"  # 0 / 20

Pinned host buffers are not affected. A host thread that keeps the GPU
busy makes the loss rare.

## First seen through MFEM

`lost-kernel.cpp` (`make lost-kernel`) zeroes components of a vector, as
`mfem::ParNonlinearForm::Mult` does for the imposed degrees of freedom,
then launches the first kernel of libmfem, whose code object is loaded at
that time (about 2 ms): the zeroes are lost. It was first taken for a
kernel lost at the first launch of a kernel of a shared library; any idle
time of the same length does it, and so do large allocations
(`hipMalloc` of 64 MiB, `hipHostMalloc` of 4 MiB).

    make lost-kernel
    sg render -c ./lost-kernel            # 0 / 27 threads, y[0] = 1
    sg render -c "./lost-kernel warm"     # libmfem loaded before: 27 / 27

## Effects in mfem-mgis

Wrong Newton residual at the first step (imposed components not zeroed),
and memory faults of hypre in rocsparse or rocprim when the lost writes
are those of hypre. The tests on this GPU run with the kernels and the
copies serialised (`AMD_SERIALIZE_KERNEL=3`, `AMD_SERIALIZE_COPY=3`,
`scripts/build-hip-amd.sh`). Even then, the examples rve and
ssna303_3d_hypre on the GPU stay blocked at the first launch of a kernel
of the library of behaviours (`BusyWaitSignal::WaitAcquire` under
`hipLaunchKernel`); not analysed further.
