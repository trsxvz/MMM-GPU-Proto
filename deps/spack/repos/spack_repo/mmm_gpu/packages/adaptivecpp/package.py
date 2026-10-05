# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import json
import os

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Adaptivecpp(CMakePackage):
    """AdaptiveCpp (formerly hipSYCL): SYCL and C++ standard parallelism
    for CPUs and GPUs, with a generic single pass compiler (SSCP).

    The `hipsycl` package of Spack stops at 24.10, which does not compile
    against CUDA 13 (`cudaDeviceProp::clockRate` was removed)."""

    homepage = "https://adaptivecpp.github.io"
    git = "https://github.com/AdaptiveCpp/AdaptiveCpp.git"

    license("BSD-2-Clause")

    provides("sycl")

    version("25.10.0", commit="9f842c701a599107cc6d117d3539f971036363a1", submodules=True)

    variant("cuda", default=False, description="Enable the CUDA backend")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.10:", type="build")
    depends_on("python@3:")
    # the released versions reject LLVM releases newer than the tested ones
    depends_on("llvm@15:20 +clang", when="@25.10.0")
    depends_on("cuda", when="+cuda")

    def cmake_args(self):
        llvm = self.spec["llvm"]
        llvm_cmake_dirs = find(llvm.prefix, "LLVMExports.cmake")
        if len(llvm_cmake_dirs) != 1:
            raise InstallError(f"LLVMExports.cmake not found in {llvm.prefix}")
        args = [
            self.define("ACPP_COMPILER_FEATURE_PROFILE", "full"),
            self.define("WITH_CUDA_BACKEND", self.spec.satisfies("+cuda")),
            self.define("WITH_ROCM_BACKEND", False),
            self.define("WITH_OPENCL_BACKEND", False),
            self.define("WITH_LEVEL_ZERO_BACKEND", False),
            self.define("LLVM_DIR", os.path.dirname(llvm_cmake_dirs[0])),
            self.define("CLANG_EXECUTABLE_PATH", llvm.prefix.bin.join("clang++")),
        ]
        if self.spec.satisfies("+cuda"):
            args.append(self.define("CUDAToolkit_ROOT", self.spec["cuda"].prefix))
            args.append(self.define("CUDA_TOOLKIT_ROOT_DIR", self.spec["cuda"].prefix))
        return args

    @run_after("install")
    def use_the_compilers_of_llvm(self):
        """The configuration of the compiler driver refers to the compiler
        used to build AdaptiveCpp, i.e. a wrapper of Spack: the clang of
        the LLVM dependency is used instead."""
        files = find(self.prefix, "acpp-core.json")
        if len(files) != 1:
            raise InstallError("acpp-core.json not found")
        with open(files[0]) as f:
            config = json.load(f)
        config["default-cpu-cxx"] = self.spec["llvm"].prefix.bin.join("clang++")
        with open(files[0], "w") as f:
            json.dump(config, f, indent=2)
