#!/usr/bin/env bash
# Launch the full pipeline detached from the terminal (survives closing it).
# Log: ../logs/full_<timestamp>.log; PID: ../logs/full.pid
# Stop it with: kill -TERM "$(cat ../logs/full.pid)"  (Snakemake resumes on the next run)
set -euo pipefail
cd "$(dirname "$0")/.."
source ../miniforge3/bin/activate pseudodefense
log="../logs/full_$(date +%Y%m%d_%H%M).log"
setsid nohup snakemake --cores 12 --rerun-incomplete > "$log" 2>&1 < /dev/null &
echo $! > ../logs/full.pid
echo "started PID $(cat ../logs/full.pid), log $log"
