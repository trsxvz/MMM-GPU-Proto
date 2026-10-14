# Kernel lost at the first launch of a kernel of a shared library

ROCm 10.1 (HIP 7.16, Python packages of TheRock), integrated GPU gfx1036
(Ryzen 7 9800X3D), MFEM 4.10 compiled with HIP (scripts/build-hip-amd.sh).

A kernel still pending on the null stream is never executed when the next
kernel is the first kernel of libmfem launched by the program, whose code
is then loaded (deferred loading, about 2 ms). The program zeroes
components of a vector, as `mfem::ParNonlinearForm::Mult` does for the
imposed degrees of freedom, then launches `other = 3.`:

    make
    sg render -c ./lost-kernel            # 0 / 27 threads, y[0] = 1
    sg render -c "./lost-kernel warm"     # libmfem warmed up: 27 / 27
    sg render -c "AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3 ./lost-kernel"

No loss with a first kernel of the executable itself. Immediate loading
(`HIP_ENABLE_DEFERRED_LOADING=0`) fails: `libhipblaslt.so.1` has no code
for gfx1036.

Effects in mfem-mgis: wrong Newton residual at the first step (imposed
components not zeroed, divergence), and memory faults of hypre in
rocsparse or rocprim when copies or fills of hypre are the ones lost.

With serialisation, some cases no longer give a wrong result but stay
blocked: in the examples rve and ssna303_3d_hypre on the GPU, the first
launch of a kernel of the library of behaviours waits forever for a
signal (`BusyWaitSignal::WaitAcquire` under `hipLaunchKernel`, stack
taken under gdb).
