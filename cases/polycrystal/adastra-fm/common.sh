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
# energy and time of the steps of this job (one step per computation)
energy() {
  sacct -j "$SLURM_JOB_ID" -P \
    --format=JobID,JobName%30,NNodes,NTasks,Elapsed,ConsumedEnergyRaw,State \
    > "$1" 2>&1
}
