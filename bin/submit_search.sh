#!/usr/bin/env bash
#SBATCH --job-name=protein_search
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/protein_search_%x_%A_%a.log
#
# SLURM launcher for protein_search.sh on large database sets (e.g. db/fungi_BFD
# has thousands of proteomes -> do not run on the head node).
#
# Splits the databases into batches and submits them as a job ARRAY (one task per
# batch, run in parallel), then submits a COLLECT job that runs after the whole
# array succeeds and builds the per-query .hits.fa files.
#
# Usage (run on the head node; forwards all protein_search.sh flags):
#   bin/submit_search.sh [-b batch_size] [-J max_concurrent] <protein_search.sh flags>
#
#   -b  databases per array task        (default: 500)
#   -J  max array tasks running at once (default: 16)
#
# Examples:
#   bin/submit_search.sh -g fungi_BFD -q fungi                 # 500/batch, phmmer
#   bin/submit_search.sh -b 250 -J 24 -g fungi_BFD -q fungi -m ssearch
#
# Resume: searches skip databases whose raw table already exists, so if some
# array tasks fail you can simply re-submit the same command to fill the gaps.
set -euo pipefail

PHASE="${SEARCH_PHASE:-launch}"

# In the launch phase we run from the real script path, so BASH_SOURCE is valid.
# Under SLURM the array/collect tasks execute a COPY of this script in the node's
# spool dir (/var/spool/slurmd/...), so BASH_SOURCE no longer points into the repo
# -- use the project paths the launcher exported instead.
if [[ "$PHASE" == "launch" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
else
  SCRIPT_DIR="$PROJECT_SCRIPT_DIR"
  ROOT="$PROJECT_ROOT"
fi
cd "$ROOT"
SEARCH="$SCRIPT_DIR/protein_search.sh"

# ---- worker: search one batch (a slice of the database manifest) -----------
if [[ "$PHASE" == "worker" ]]; then
  start=$(( SLURM_ARRAY_TASK_ID * BATCH + 1 ))
  # shellcheck disable=SC2086  # PASS is an intentional flag string
  exec "$SEARCH" $PASS -L "$MANIFEST" -r "${start}:${BATCH}" -s -t "${SLURM_CPUS_PER_TASK:-8}"
fi

# ---- collect: build .hits.fa from all raw tables ---------------------------
if [[ "$PHASE" == "collect" ]]; then
  # shellcheck disable=SC2086
  exec "$SEARCH" $PASS -c
fi

# ---- launcher (runs on the head node) --------------------------------------
# Collect everything that isn't -b/-J into PASS: the protein_search.sh flags to
# forward. On the head node PASS is an array; SLURM --export can only carry
# strings, so below we flatten it with "${PASS[*]}" and the worker/collect phases
# re-split the unquoted $PASS back into separate flags.
BATCH=500
MAXJOBS=16
PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) BATCH="$2"; shift 2 ;;
    -J) MAXJOBS="$2"; shift 2 ;;
    *)  PASS+=("$1"); shift ;;
  esac
done

mkdir -p logs
MANIFEST="$ROOT/logs/dblist.$$.txt"
"$SEARCH" "${PASS[@]}" -p > "$MANIFEST"
N=$(wc -l < "$MANIFEST")
[[ "$N" -gt 0 ]] || { echo "ERROR: no databases found" >&2; exit 1; }
NTASKS=$(( (N + BATCH - 1) / BATCH ))
echo "==> $N databases -> $NTASKS array tasks of up to $BATCH (max $MAXJOBS concurrent)"
echo "==> manifest: $MANIFEST"

# search array
AJ=$(sbatch --parsable --job-name=psearch \
     --array="0-$((NTASKS-1))%${MAXJOBS}" \
     --export="ALL,SEARCH_PHASE=worker,BATCH=$BATCH,MANIFEST=$MANIFEST,PASS=${PASS[*]},PROJECT_SCRIPT_DIR=$SCRIPT_DIR,PROJECT_ROOT=$ROOT" \
     "$SCRIPT_DIR/submit_search.sh")
echo "==> search array job: $AJ"

# collect after the whole array completes successfully
CJ=$(sbatch --parsable --job-name=pcollect \
     --dependency="afterok:$AJ" --cpus-per-task=2 --array=0 \
     --export="ALL,SEARCH_PHASE=collect,PASS=${PASS[*]},PROJECT_SCRIPT_DIR=$SCRIPT_DIR,PROJECT_ROOT=$ROOT" \
     "$SCRIPT_DIR/submit_search.sh")
echo "==> collect job (after array): $CJ"
echo "    monitor: squeue -j $AJ,$CJ ; logs in logs/"
