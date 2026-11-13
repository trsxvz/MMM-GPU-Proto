#!/bin/bash
# Libraries of the variants of the behaviour of the polycrystal
# (polycrystal/law-variants/Mono_UO2_CosH_Jaco3-<variant>.mfront, copied in
# <scratch>/laws/<variant>/), built as the ones of the prototype: for the
# GPU by the generic-parallel interface (hip, as cases/polycrystal), for the
# CPU by the generic interface (as mm-opera-hpc/polycrystal), on a node of
# the MI300A partition. FM_LAWS: variants (default "A B").
#SBATCH --job-name=fm-build-laws
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus-per-node=4
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
mfront=$D/MMM-GPU-Proto/install/tfel/bin/mfront
log "start, variants ${FM_LAWS:-A B}"
for v in ${FM_LAWS:-A B}; do
  l=$W/laws/$v
  mkdir -p "$l/gpu" "$l/cpu"
  (cd "$l/gpu" && "$mfront" --obuild --interface=generic-parallel \
    --configuration-file=config-hip-hipcc.json ../Mono_UO2_CosH_Jaco3.mfront) \
    > "$l/gpu.log" 2>&1
  (cd "$l/cpu" && "$mfront" --obuild --interface=generic \
    ../Mono_UO2_CosH_Jaco3.mfront) > "$l/cpu.log" 2>&1
  log "$v: $(ls "$l"/gpu/src/libGenericParallelBehaviour-hip.so \
    "$l"/cpu/src/libBehaviour.so 2>&1 | tr '\n' ' ')"
done
log end
