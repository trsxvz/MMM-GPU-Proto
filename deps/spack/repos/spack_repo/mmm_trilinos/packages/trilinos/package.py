# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.trilinos.package import Trilinos as BuiltinTrilinos

from spack.package import *


class Trilinos(BuiltinTrilinos):
    """Trilinos of the builtin repository, with the sparse direct solver
    Tacho in Amesos2 (subpackage ShyLU_NodeTacho, on the host and on the
    GPU), which the package of spack does not enable: the variant shylu
    only enables ShyLU_DD (FROSch), and the optional dependencies of Amesos2
    are not enabled by default."""

    variant(
        "tacho",
        default=False,
        when="+amesos2 +shylu",
        description="Sparse direct solver Tacho in Amesos2",
    )

    def cmake_args(self):
        args = super().cmake_args()
        if self.spec.satisfies("+tacho"):
            # the tests of Tacho require GoogleTest, which TriBITS
            # requires as soon as Tacho is enabled, even without the tests
            args += [
                self.define("Trilinos_ENABLE_ShyLU_NodeTacho", True),
                self.define("Amesos2_ENABLE_ShyLU_NodeTacho", True),
                self.define("ShyLU_NodeTacho_ENABLE_TESTS", False),
                self.define("TPL_ENABLE_gtest", False),
            ]
        return args
