#!/bin/bash
# Run the Gaia XP mirror on one whole CPU node: NWORK workers, each a fixed round-robin share of
# the 3,386 CDN files (resumable through _ledger/gaia_xp.jsonl).  Started from a login node as
#   nohup setsid salloc -A <acct> -q interactive -C cpu -N 1 -t 04:00:00 \
#        srun -N 1 -n 1 -c 256 --cpu-bind=none bash slurm/gaia_xp_node.sh > _logs/gaia_xp_node.log 2>&1 &
# (~8 GB peak memory per worker: 32 workers fit a 512 GB node.)
NWORK=${NWORK:-32}
module load python
export MIRROR_STAGING=$PSCRATCH/dustline_mirror_staging${XP_TAG}
cd /global/cfs/projectdirs/newera/surveys/_tools/mirror
LOG=/global/cfs/projectdirs/newera/surveys/_logs
echo "node $(hostname), $NWORK workers, start $(date)"
for i in $(seq 0 $((NWORK - 1))); do
    python gaia_xp.py run --task-id "$i" --ntasks "$NWORK" $XP_ARGS > "$LOG/gaia_xp${XP_TAG}_w$(printf %02d "$i").log" 2>&1 &
done
wait
echo "all workers finished $(date)"
python gaia_xp.py manifest
