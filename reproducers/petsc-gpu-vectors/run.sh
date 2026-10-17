#!/bin/bash
# Build the example ex1p of PETSc of MFEM 4.10 (unmodified copy) against the
# MFEM of the prototype (CUDA, PETSc 3.26 built with CUDA) and solve the
# same problem on the CPU and on the GPU, on 1 and 2 processes: on the GPU,
# the preconditioned residual of the conjugate gradient of PETSc is zero
# from the first iteration and the solution stays zero.
# usage: reproducers/petsc-gpu-vectors/run.sh [MFEM installation]
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
M=${1:-$(grep "^MFEM_DIR" "$P/build/bench/mfem-mgis/CMakeCache.txt" |
         cut -d= -f2)}
MPIRUN=$(dirname "$(grep "^MFEM_HOST_CXX" "$M/share/mfem/config.mk" |
                    awk '{print $3}')")/mpirun
W=$P/scratch/reproducers/petsc-gpu-vectors
mkdir -p "$W" && cd "$W" || exit 1
cp "$R"/{ex1p.cpp,makefile,rc_ex1p_device,star.mesh} . &&
  make -s MFEM_DIR="$M" ex1p > build.log 2>&1 || { tail build.log; exit 1; }
export OMPI_MCA_btl=self,tcp
for args in "1 cpu" "1 cuda" "2 cpu" "2 cuda" "2 cuda -pa"; do
  set -- $args
  d=run-$(echo "$args" | tr ' ' '-')
  mkdir -p "$d" && cp star.mesh rc_ex1p_device "$d" || exit 1
  (cd "$d" && PETSC_OPTIONS="-ksp_monitor_true_residual" timeout 120 \
     "$MPIRUN" -n "$1" ../ex1p -m star.mesh -no-vis -d "$2" "${@:3}" \
     --petscopts rc_ex1p_device > out.log 2>&1)
  echo "== $args: exit $?"
  grep -m1 " 0 KSP" "$d/out.log" | cut -c1-110
  grep -E "Linear solve" "$d/out.log"
  python3 - "$d/sol.000000" <<'PY'
import sys
values = []
for line in open(sys.argv[1]).read().split("\n\n")[-1].split():
    try:
        values.append(abs(float(line)))
    except ValueError:
        pass
print(f"  largest value of the solution: {max(values):.6g}")
PY
done
