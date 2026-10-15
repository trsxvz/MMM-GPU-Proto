# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.mmm_gpu.packages.petsc.package import Petsc as GpuPetsc

from spack.package import *


class Petsc(GpuPetsc):
    """PETSc of the repository mmm_gpu, with the sparse direct solver
    PaStiX (6.4, the version supported by PETSc 3.26), which the package of
    spack does not enable. PaStiX requires hwloc and an ordering library
    (METIS here)."""

    variant("pastix", default=False, when="+mpi", description="Activates support for PaStiX")

    depends_on("pastix@6.4: +mpi", when="+pastix")
    conflicts("~hwloc", when="+pastix", msg="PaStiX requires hwloc")
    conflicts("~metis", when="+pastix", msg="PaStiX requires an ordering library (METIS)")

    def configure_options(self):
        options = super().configure_options()
        if self.spec.satisfies("+pastix"):
            options += ["--with-pastix=1", "--with-pastix-dir=%s" % self.spec["pastix"].prefix]
        else:
            options += ["--with-pastix=0"]
        return options
