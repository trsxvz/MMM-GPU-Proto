# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.amgx.package import Amgx as BuiltinAmgx

from spack.package import *


class Amgx(BuiltinAmgx):
    """AmgX, with a commit of the main branch compatible with CUDA 13.

    AmgX 2.5.0 predates the fixes for CUDA 13 (C++17 standard with CUDA 13.1
    and later, compatibility of the tuples of CUDA 13), made on the main
    branch in 2026."""

    git = "https://github.com/NVIDIA/AMGX.git"

    version("2.5.0.20260701", commit="91a8413ef267b1c32aff4014c02820e1c5897ac2")
