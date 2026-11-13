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
log() { echo "$(date '+%F %T') job $SLURM_JOB_ID ($SLURM_JOB_NAME): $*" >> "$R/chain.log"; }
# energy and time of the steps of this job (one step per computation)
energy() {
  sacct -j "$SLURM_JOB_ID" -P \
    --format=JobID,JobName%30,NNodes,NTasks,Elapsed,ConsumedEnergyRaw,State \
    > "$1" 2>&1
}
