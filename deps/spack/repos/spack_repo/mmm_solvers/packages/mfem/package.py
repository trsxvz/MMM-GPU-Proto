# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.mmm_gpu.packages.mfem.package import Mfem as GpuMfem

from spack.package import *


class Mfem(GpuMfem):
    """MFEM of the repository mmm_gpu, with the sparse direct solvers of
    oneMKL (PARDISO, sequential, and the Cluster Sparse Solver, parallel)
    and MAGMA (batched dense linear algebra on the GPU), which the package
    of spack does not enable.

    The Cluster Sparse Solver uses the BLACS of oneMKL for the MPI library
    of the environment (OpenMPI or an MPI of the family of MPICH). oneMKL is
    given by the package intel-oneapi-mkl-solvers of mmm_solvers, which does
    not provide BLAS and LAPACK (see this package)."""

    variant("mkl-pardiso", default=False, description="Enable the PARDISO solver of oneMKL")
    variant(
        "mkl-cpardiso",
        default=False,
        when="+mpi",
        description="Enable the Cluster Sparse Solver (parallel PARDISO) of oneMKL",
    )
    variant("magma", default=False, description="Enable MAGMA (batched dense linear algebra)")

    depends_on("intel-oneapi-mkl-solvers", when="+mkl-pardiso")
    # with the BLACS of oneMKL
    depends_on("intel-oneapi-mkl-solvers", when="+mkl-cpardiso")
    depends_on("magma", when="+magma")
    depends_on("magma +cuda", when="+magma +cuda")
    for arch in CudaPackage.cuda_arch_values:
        depends_on(f"magma cuda_arch={arch}", when=f"+magma +cuda cuda_arch={arch}")

    def get_make_config_options(self, spec, prefix):
        options = super().get_make_config_options(spec, prefix)

        def yes_no(v):
            return "YES" if spec.satisfies(v) else "NO"

        options += [
            "MFEM_USE_MKL_PARDISO=%s" % yes_no("+mkl-pardiso"),
            "MFEM_USE_MKL_CPARDISO=%s" % yes_no("+mkl-cpardiso"),
            "MFEM_USE_MAGMA=%s" % yes_no("+magma"),
        ]
        if spec.satisfies("+mkl-pardiso") or spec.satisfies("+mkl-cpardiso"):
            # the directory of oneMKL, with include and lib
            mkl_dir = os.path.dirname(spec["intel-oneapi-mkl-solvers"].headers.directories[0])
            mkl_lib = os.path.join(mkl_dir, "lib")
            core = ["mkl_intel_lp64", "mkl_sequential", "mkl_core"]
            if spec.satisfies("+mkl-pardiso"):
                options += [
                    "MKL_PARDISO_DIR=%s" % mkl_dir,
                    "MKL_PARDISO_OPT=-I%s" % os.path.join(mkl_dir, "include"),
                    "MKL_PARDISO_LIB=%s" % self.ld_flags_from_dirs([mkl_lib], core),
                ]
            if spec.satisfies("+mkl-cpardiso"):
                blacs = (
                    "mkl_blacs_openmpi_lp64"
                    if spec.satisfies("^[virtuals=mpi] openmpi")
                    else "mkl_blacs_intelmpi_lp64"
                )
                options += [
                    "MKL_CPARDISO_DIR=%s" % mkl_dir,
                    "MKL_CPARDISO_OPT=-I%s" % os.path.join(mkl_dir, "include"),
                    "MKL_CPARDISO_LIB=%s" % self.ld_flags_from_dirs([mkl_lib], [blacs] + core),
                ]
        if spec.satisfies("+magma"):
            magma = spec["magma"]
            options += [
                "MAGMA_DIR=%s" % magma.prefix,
                "MAGMA_OPT=-I%s" % magma.prefix.include,
                "MAGMA_LIB=%s"
                % self.ld_flags_from_library_list(magma.libs + spec["lapack"].libs),
            ]
        return options
