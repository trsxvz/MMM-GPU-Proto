#!/bin/bash
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
# Archives in parallel the time steps (CycleNNNNNN) of ParaView exports, one
# archive per step (archive-cycle.sh), a process per step spread over the
# nodes of the allocation.
# usage: sbatch -A <account> --nodes=<n> --time=<t> \
#          --export=ALL,ARCHIVE_DIRS="<export directory> ..." archive-cycles.slurm
set -u
here=${ARCHIVE_SCRIPTS:-$SLURM_SUBMIT_DIR/fm-chain}
# the list is read by the processes of every node: on Lustre, not in /tmp
list=$(mktemp "$SLURM_SUBMIT_DIR/.cycles-$SLURM_JOB_ID-XXXX")
find ${ARCHIVE_DIRS:?} -maxdepth 1 -type d -name 'Cycle*' | sort > "$list"
n=$(wc -l < "$list")
[ "$n" -gt 0 ] || exit 0
start=$(date +%s)
srun --ntasks="$n" --nodes="$(( n < SLURM_JOB_NUM_NODES ? n : SLURM_JOB_NUM_NODES ))" \
  --cpus-per-task=1 --distribution=cyclic bash -c \
  'bash "$0" "$(sed -n "$((SLURM_PROCID + 1))p" "$1")"' "$here/archive-cycle.sh" "$list"
echo "$n steps archived in $(( $(date +%s) - start )) s:" \
  "$(find ${ARCHIVE_DIRS} -maxdepth 1 -type d -name 'Cycle*' | wc -l) left unarchived"
rm -f "$list"
