# Sourced by the jobs of the chain (fm-chain/*.slurm): environment of the
# prototype on the MI300A partition of Adastra, results directory R.
set -u
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ftp_proxy FTP_PROXY
D=/lus/work/RES1/genden15/tchenaille/mmm
export PROTO_DEPLOY=$D PROTO_MACHINE=adastra-mi300
# prototype: the deployed one, or a version of its scripts (FM_PROTO, see
# deploy-overlay.sh)
P=${FM_PROTO:-$D/MMM-GPU-Proto}
# results and meshes on the scratch: the quota of files of the project on
# /lus/work is reached (2026-11-13)
W=/lus/scratch/RES1/genden15/tchenaille
R=$W/fm/${FM_RUN:?FM_RUN required}
MESHES=$W/meshes/series
mkdir -p "$R"
# journal of the jobs: chain.log of the results directory and, if given, a
# second journal FM_CHAIN_LOG (for instance the one of the campaign)
log() {
  echo "$(date '+%F %T') job $SLURM_JOB_ID ($SLURM_JOB_NAME): $*" |
    tee -a ${FM_CHAIN_LOG:+"$FM_CHAIN_LOG"} >> "$R/chain.log"
}
# chain of jobs: the job stops at once unless the files FM_REQUIRE exist
# (files written by the previous jobs of the chain when their checks pass)
for _required in ${FM_REQUIRE:-}; do
  if [ ! -e "$_required" ]; then
    log "required file $_required missing: job stopped"
    exit 1
  fi
done
unset _required
# usage: load_build_environment
# environment of the deployed prototype for the jobs which build or run
# tests (compilers, libraries, mpirun of the machine profile). env.sh
# redefines P, the deployed prototype: the scripts of the computations
# remain the ones of FM_PROTO
load_build_environment() {
  local scripts=$P
  set +u
  if ! source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null; then
    log "environment of the prototype not loaded"
    exit 1
  fi
  set -u
  P=$scripts
  export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
  unset OMP_PROC_BIND OMP_PLACES
}
# usage: behaviour_targets <build directory of mfem-mgis>
# targets of the behaviours of the tests, which the executables of the
# tests do not build
behaviour_targets() {
  echo BehaviourTest
  cmake --build "$1" --target help |
    sed -n 's/^\.\.\. \(GenericParallelBehaviour-[a-z]*\)$/\1/p'
}
# usage: test_settings mode
# variables of a run of the tests of mfem-mgis by ctest, one per line:
# - host: on the host, one thread per process;
# - host-kernels: on the host, all the computations forced to the kernels;
# - gpu-kernels, gpu-default: on a GPU, with the settings of
#   scripts/mfem-mgis-gpu-suite.sh (device of the machine profile, one GiB
#   for the pool of hypre, stiffness matrices added in a fixed order), the
#   computations being forced to the kernels or chosen by each test
test_settings() {
  local kernels="MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels"
  local gpu="MFEM_MGIS_DEVICE=${PROTO_GPU_DEVICE:-cuda}
MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1
MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered"
  case $1 in
    host) echo OMP_NUM_THREADS=1 ;;
    host-kernels) printf '%s\n' OMP_NUM_THREADS=1 "$kernels" ;;
    gpu-kernels) printf '%s\n' "$gpu" "$kernels" ;;
    gpu-default) echo "$gpu" ;;
    *) return 1 ;;
  esac
}
# energy and time of the steps of this job (one step per computation)
energy() {
  sacct -j "$SLURM_JOB_ID" -P \
    --format=JobID,JobName%30,NNodes,NTasks,Elapsed,ConsumedEnergyRaw,State \
    > "$1" 2>&1
}
