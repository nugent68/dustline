#!/bin/bash
# Gaia DR3 XP continuous mirror: download (md5-checked), convert to FITS, verify, delete raw.
# Resumable: resubmit the same array to pick up anything that failed (see _ledger/gaia_xp.jsonl).
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH -c 4
#SBATCH --mem=14G
#SBATCH -t 04:00:00
#SBATCH -J gaia_xp
#SBATCH --array=0-63%32
#SBATCH -o /global/cfs/projectdirs/newera/surveys/_logs/gaia_xp_%A_%a.log

module load python
export MIRROR_STAGING=$PSCRATCH/dustline_mirror_staging
cd /global/cfs/projectdirs/newera/surveys/_tools/mirror
python gaia_xp.py run --task-id "$SLURM_ARRAY_TASK_ID" --ntasks 64
