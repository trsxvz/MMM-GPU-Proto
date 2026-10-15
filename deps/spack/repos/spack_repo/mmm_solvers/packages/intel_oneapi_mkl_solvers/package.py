# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.oneapi import (
    IntelOneApiLibraryPackage,
    IntelOneApiPackage,
)

from spack.package import *


@IntelOneApiPackage.update_description
class IntelOneapiMklSolvers(IntelOneApiLibraryPackage):
    """oneMKL installed for its sparse direct solvers (PARDISO and the
    Cluster Sparse Solver), with the installer of the package
    intel-oneapi-mkl, but without providing BLAS, LAPACK and ScaLAPACK: in a
    unified environment of spack, a package which can provide a virtual
    package present in the graph is its provider for the whole graph, and
    oneMKL would replace the BLAS and LAPACK of the environment (OpenBLAS)
    for all the packages."""

    homepage = (
        "https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html"
    )

    version(
        "2025.3.1",
        url="https://registrationcenter-download.intel.com/akdlm/IRC_NAS/6a17080f-f0de-41b9-b587-52f92512c59a/intel-onemkl-2025.3.1.11_offline.sh",
        sha256="89753ce0be82d31669172c08c6b6b863f2e25558d775496349473b3299240a01",
        expand=False,
    )

    @property
    def v2_layout_versions(self):
        return "@2024:"

    @property
    def component_dir(self):
        return "mkl"

    @property
    def headers(self):
        headers = find_all_headers(self.component_prefix.include)
        headers.directories = [self.component_prefix.include]
        return headers

    @property
    def libs(self):
        # the libraries of the sparse direct solvers, sequential, LP64
        return find_libraries(
            ["libmkl_intel_lp64", "libmkl_sequential", "libmkl_core"],
            self.component_prefix.lib,
            shared=True,
        )
