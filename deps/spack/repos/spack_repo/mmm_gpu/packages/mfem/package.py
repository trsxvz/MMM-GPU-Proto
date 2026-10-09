# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import filecmp
import os
import shutil

from spack_repo.builtin.packages.mfem.package import Mfem as BuiltinMfem

from spack.package import *


class Mfem(BuiltinMfem):
    """MFEM, whose partial builds are resumed where they stopped.

    The objects of MFEM depend on `config/config.mk`, written again by
    `make config` at each build: resuming a partial build (spack install
    --dont-restage) would compile everything again (several hours with CUDA).
    The previous configuration files keep their date when their content does
    not change."""

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
