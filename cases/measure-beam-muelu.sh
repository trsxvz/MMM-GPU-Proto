#!/bin/bash
# Newton solver on the beam (cases/beam) with the multigrid of MueLu
# (conjugate gradient of MFEM, hierarchy reused from one tangent operator to
# the next) against BoomerAMG of hypre, on the GPU (one process) and on the
# CPU (8 processes). The cases must be built against an installation of
# mfem-mgis built with Trilinos.
#
# usage: measure-beam-muelu.sh [cases build] [results directory]
#   cases build: by default build/cases-trilinos;
#   results directory: by default results/beam-muelu-<date>.
# Environment: REFINEMENT (4 by default: 111,843 unknowns), TARGETS ("gpu
# cpu" by default), FILTER (regular expression selecting the linear
# solvers), OPTIONS (other options of the beam, for instance `-pu 4`: the
# preconditioner of cg-amg and cg-muelu built every 4 Newton iterations)
# and LABEL (suffix of the names of the runs).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
C=$(cd "${1:-$P/build/cases-trilinos}" && pwd)
out=${2:-$P/results/beam-muelu-$(date +%Y-%m-%d-%H%M)}
B=$C/beam
refinement=${REFINEMENT:-4}
targets=${TARGETS:-gpu cpu}
filter=${FILTER:-.}
options=${OPTIONS:-}
label=${LABEL:-}
if [ ! -x "$B/beam" ]; then
  echo "no beam case in $B" >&2
  exit 1
fi
mkdir -p "$out"
MPI=$(dirname "$(dirname "$(/usr/bin/grep -m1 '^MPI_CXX_COMPILER:' "$C/CMakeCache.txt" | sed 's/.*=//')")")
CUDA=${PROTO_CUDA:-$HOME/Bureau/spack-mmm-externals/cuda-13.0.2}
T=$P/install/h-interface/tfel
M=$P/install/h-interface/mgis
export PATH=$CUDA/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=$T/lib:$M/lib:$CUDA/lib64
export OMPI_MCA_btl=self,vader
export OMP_NUM_THREADS=1
unset OMP_PROC_BIND OMP_PLACES
# behaviours integrated by kernels
export MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
export MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
export MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
export MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
# context and summary appended, the script being run for several
# refinements in the same directory
{
  echo "date: $(date '+%F %T')"
  echo "cases: ${C#"$P"/}"
  echo "refinement: $refinement"
  echo "mfem-mgis: $(/usr/bin/grep -m1 '^MFEMMGIS_DIR:' "$C/CMakeCache.txt" | sed "s|.*=||; s|$P/||")"
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null
  LC_ALL=C lscpu | /usr/bin/grep "Model name"
} >> "$out/context.txt"
for target in $targets; do
  if [ "$target" = gpu ]; then
    launcher=(taskset -c 0,1 mpirun -np 1 --bind-to none -x LD_LIBRARY_PATH
      -x OMPI_MCA_btl -x OMP_NUM_THREADS -x MFEM_MGIS_BEHAVIOUR_INTEGRATION
      -x MFEM_MGIS_GRADIENTS_COMPUTATION -x MFEM_MGIS_INNER_FORCES_COMPUTATION
      -x MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION)
    library=$B/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
    device=(--device cuda)
    solvers="pcg-amg cg-amg cg-muelu"
  else
    launcher=(taskset -c 0-7 mpirun -np 8 --bind-to core -x LD_LIBRARY_PATH
      -x OMPI_MCA_btl -x OMP_NUM_THREADS -x MFEM_MGIS_BEHAVIOUR_INTEGRATION
      -x MFEM_MGIS_GRADIENTS_COMPUTATION -x MFEM_MGIS_INNER_FORCES_COMPUTATION
      -x MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION)
    library=$B/src/libBehaviour.so
    device=()
    solvers="pcg-amg cg-muelu-sgs cg-muelu"
  fi
  for ls in $solvers; do
    echo "$ls" | /usr/bin/grep -Eq "$filter" || continue
    log=$out/beam-r$refinement-$target-$ls$label.txt
    (cd "$B" && timeout 1800 "${launcher[@]}" ./beam -r "$refinement" -n 10 \
      -mi 3000 -ls "$ls" -lt 1e-4 -l "$library" --profiling "${device[@]}" \
      $options) \
      < /dev/null > "$log" 2>&1
    echo "r$refinement $target $ls$label: status $? $(/usr/bin/grep -E '^dofs|^time' "$log" | tr '\n' ' ')" | tee -a "$out/summary.txt"
  done
done
