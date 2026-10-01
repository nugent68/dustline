#!/bin/bash
# One whole CPU node: AllWISE column extract (reads the 1.5 TB cosmo copy), then - once the Dataverse
# download has all 100 batches - the Zucker+25 stellar_inference conversion (~200 GB memory).
#   nohup setsid salloc -A <acct> -q interactive -C cpu -N 1 -t 04:00:00 \
#        srun -N 1 -n 1 -c 256 --cpu-bind=none bash slurm/node_allwise_zucker.sh > _logs/node_allwise_zucker.log 2>&1 &
module load python
S=/global/cfs/projectdirs/newera/surveys
export DUSTLINE_SURVEYS=$S
cd $S/_tools/mirror
echo "start $(date) on $(hostname)"
python extract.py allwise --nproc 24 && echo "allwise done $(date)"
until [ "$(grep -c '^ok' $S/_logs/zucker25_download.log)" -ge 100 ]; do
    grep -q FAILED $S/_logs/zucker25_download.log && { echo "zucker download has FAILED files"; exit 1; }
    sleep 60
done
python zucker25.py convert && python zucker25.py manifest && echo "zucker done $(date)"
