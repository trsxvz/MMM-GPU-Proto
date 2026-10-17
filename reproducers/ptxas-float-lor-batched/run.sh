#!/bin/bash
# Memory of ptxas on fem/lor/lor_batched.cpp of MFEM 4.10 built in single
# precision (MFEM_USE_SINGLE=YES) with CUDA 13.0.2 for sm_120. The command
# of the build of MFEM (nvcc-command.txt, recorded from the spack build of
# the environment deps/spack/mmm-trilinos-float) is run in the source tree
# of this build, given as first argument, with the target given as second
# argument (sm_120: PTX and machine code, ptxas run; compute_120: PTX only)
# and the ptxas options given as third argument. The largest resident size
# among the processes of the compilation is printed by GNU time.
# usage: run.sh mfem-source-tree [sm_120|compute_120] [ptxas-options]
set -eu
here=$(cd "$(dirname "$0")" && pwd)
src=${1:?usage: run.sh mfem-source-tree [sm_120|compute_120] [ptxas-options]}
target=${2:-sm_120}
ptxas=${3:-}
command=$(sed "s| -arch=sm_120 | -arch=$target ${ptxas:+-Xptxas $ptxas }|; s| -o fem/lor/lor_batched.o| -o /dev/null|" "$here/nvcc-command.txt")
cd "$src"
/usr/bin/time -f "maximum resident size: %M kB, wall time: %e s" bash -c "$command"
