# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.openmpi.package import Openmpi as BuiltinOpenmpi

from spack.package import *


class Openmpi(BuiltinOpenmpi):
    """OpenMPI, whose internal hwloc is built without its OpenCL plugin.

    With CUDA, the internal hwloc of OpenMPI 5 finds the OpenCL library of
    the CUDA toolkit, enables its OpenCL plugin and fails to link it
    (`cannot determine absolute directory name of 'no/opencl/lib/'`). The
    hwloc of spack cannot be used instead: hwloc with CUDA, shared with LLVM
    and TBB, would rebuild them."""

    def configure_args(self):
        args = super().configure_args()
        if self.spec.satisfies("+internal-hwloc"):
            args.append("--disable-opencl")
        return args
