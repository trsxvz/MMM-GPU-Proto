# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.amgx.package import Amgx as BuiltinAmgx

from spack.package import *


class Amgx(BuiltinAmgx):
    """AmgX, with a commit of the main branch compatible with CUDA 13.

    AmgX 2.5.0 predates the fixes for CUDA 13 (C++17 standard with CUDA 13.1
    and later, compatibility of the tuples of CUDA 13), made on the main
    branch in 2026.

    With CUDA 13, the device code of src/amg.cu lacks the kernel of
    thrust::sequence on the device vectors of integers, whose host code is
    instantiated there by the inline functions of the matrices: the copy of
    cub::DeviceFor::Bulk of this file launches a null kernel, and, kept by
    the linker, makes the setup of the classical multigrid fail with
    "invalid device function". The defect does not depend on the option
    -static-global-template-stub of nvcc (the null launch remains with
    =false). The patch uses thrust::sequence in a function of this file,
    whose device code then contains the kernel (reproducer:
    reproducers/amgx-classical-cuda13 of the prototype)."""

    git = "https://github.com/NVIDIA/AMGX.git"

    version("2.5.0.20260701", commit="91a8413ef267b1c32aff4014c02820e1c5897ac2")

    patch("amgx-cuda13-sequence-kernel.patch", when="@2.5.0.20260701 ^cuda@13:")
