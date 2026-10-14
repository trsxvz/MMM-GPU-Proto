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
    The previous configuration files keep their date when the lines which
    matter for the compilation of the objects do not change.

    With CUDA 13, `-lnvToolsExt` is removed from the link flags: the default
    link flags of AmgX in MFEM name this library, removed from CUDA 13 (NVTX 3
    is header-only, and AmgX only uses its headers).

    With CUDA 13 (sm_120), cicc does not finish compiling
    `fem/integ/bilininteg_convection_ea.cpp` (more than 2 h 30 and 17 GB,
    then killed for lack of memory): its 3D kernels are specialized up to 8
    degrees of freedom and 9 quadrature points per direction. The patch keeps
    the specializations up to the order 3 and leaves the higher orders to the
    generic kernel of the same file.

    The variant `cuda_ptx_only` embeds the PTX code only, without SASS:
    ptxas is not run at compilation, the driver compiles the PTX code at
    the first run. It is meant for the compilation checks of the
    continuous integration: with CUDA 13 and sm_120, ptxas needs more than
    12 GB for `fem/lor/lor_batched.cpp`, more than the memory of the
    machines of GitHub."""

    patch("mfem-4.10-convection-ea-cuda13.patch", when="@4.10.0 +cuda ^cuda@13:")

    variant(
        "cuda_ptx_only",
        default=False,
        when="+cuda",
        description="Embed the PTX code only (no SASS, ptxas not run at "
        "compilation): for compilation checks",
    )

    _configuration_files = ["config/config.mk", "config/_config.hpp"]

    @staticmethod
    def _compilation_lines(f):
        """lines of a configuration file which matter for the compilation of
        the objects: the libraries given to the linker (`MFEM_EXT_LIBS`)
        do not"""
        with open(f) as fh:
            return [line for line in fh if not line.startswith("MFEM_EXT_LIBS")]

    def configure(self, spec, prefix):
        with working_dir(self.stage.source_path):
            for f in self._configuration_files:
                if os.path.exists(f):
                    shutil.copy2(f, f + ".mmm-previous")
        super().configure(spec, prefix)
        with working_dir(self.stage.source_path):
            if spec.satisfies("+cuda ^cuda@13:"):
                # the default link flags of AmgX (config/defaults.mk) name
                # libnvToolsExt, removed from CUDA 13 (NVTX 3 is header-only)
                filter_file(r" -lnvToolsExt\b", "", "config/config.mk")
            if spec.satisfies("+cuda_ptx_only"):
                # virtual architecture only: the PTX code is embedded and
                # ptxas is not run
                filter_file(
                    r"-arch=sm_(\d+)",
                    r"-arch=compute_\1 -code=compute_\1",
                    "config/config.mk",
                )
            for f in self._configuration_files:
                previous = f + ".mmm-previous"
                if not os.path.exists(previous):
                    continue
                if os.path.exists(f) and (
                    self._compilation_lines(f) == self._compilation_lines(previous)
                ):
                    s = os.stat(previous)
                    os.utime(f, (s.st_atime, s.st_mtime))
                os.remove(previous)
