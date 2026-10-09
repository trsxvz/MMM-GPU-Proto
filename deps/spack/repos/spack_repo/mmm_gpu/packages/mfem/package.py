# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import filecmp
import os
import shutil

from spack_repo.builtin.packages.mfem.package import Mfem as BuiltinMfem

from spack.package import *


class Mfem(BuiltinMfem):
    """MFEM, whose partial builds are resumed where they stopped, and whose
    kernels of high order of the element assembly of the convection
    integrator in 3D are not specialized with CUDA.

    The objects of MFEM depend on `config/config.mk`, written again by
    `make config` at each build: resuming a partial build (spack install
    --dont-restage) would compile everything again (several hours with CUDA).
    The previous configuration files keep their date when their content does
    not change.

    With CUDA 13 (sm_120), cicc does not finish compiling
    `fem/integ/bilininteg_convection_ea.cpp` (more than 2 h 30 and 17 GB,
    then killed for lack of memory): its 3D kernels are specialized up to 8
    degrees of freedom and 9 quadrature points per direction. The patch keeps
    the specializations up to the order 3 and leaves the higher orders to the
    generic kernel of the same file."""

    patch("mfem-4.10-convection-ea-cuda13.patch", when="@4.10.0 +cuda")

    _configuration_files = ["config/config.mk", "config/_config.hpp"]

    def configure(self, spec, prefix):
        with working_dir(self.stage.source_path):
            for f in self._configuration_files:
                if os.path.exists(f):
                    shutil.copy2(f, f + ".mmm-previous")
        super().configure(spec, prefix)
        with working_dir(self.stage.source_path):
            for f in self._configuration_files:
                previous = f + ".mmm-previous"
                if not os.path.exists(previous):
                    continue
                if os.path.exists(f) and filecmp.cmp(f, previous, shallow=False):
                    s = os.stat(previous)
                    os.utime(f, (s.st_atime, s.st_mtime))
                os.remove(previous)
