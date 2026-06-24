#!/usr/bin/env bash
#
# protein_search.sh - search query proteins against grouped sequence databases
#                     and collect high-scoring hits into per-query multi-FASTA files.
#
# Layout assumed:
#   <group>/*.fa[a|sta]        query proteins for a major group (e.g. fungi/, insect/)
#   db/<group>/*.fasta         subject databases for that major group
#   results/<group>/           output (created here)
#
# For each database in db/<group>/, every query in <group>/ is searched once.
# Hits with E-value < threshold are fetched from the source database and written
# to results/<group>/<query_id>.hits.fa, one multi-FASTA per query protein.
#
# Usage:
#   bin/protein_search.sh [-g group] [-m method] [-e evalue] [-t threads]
#
#   -g  major group / query+db folder name   (default: fungi)
#   -m  search method: ssearch | phmmer      (default: phmmer)
#   -e  E-value threshold (hits must be <)    (default: 1e-20)
#   -t  CPU threads per search               (default: 8)
#
# Examples:
#   bin/protein_search.sh                          # fungi, phmmer, 1e-20
#   bin/protein_search.sh -g fungi -m ssearch      # Smith-Waterman via SSEARCH
#   bin/protein_search.sh -g insect -e 1e-30 -t 16
#
set -euo pipefail

# ---- defaults --------------------------------------------------------------
GROUP="fungi"
METHOD="phmmer"
EVALUE="1e-20"
THREADS="8"

while getopts "g:m:e:t:h" opt; do
  case "$opt" in
    g) GROUP="$OPTARG" ;;
    m) METHOD="$OPTARG" ;;
    e) EVALUE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Run '$0 -h' for usage." >&2; exit 1 ;;
  esac
done

# ---- resolve project root (parent of bin/) --------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

QUERY_DIR="$ROOT/$GROUP"
DB_DIR="$ROOT/db/$GROUP"
OUT_DIR="$ROOT/results/$GROUP"
RAW_DIR="$OUT_DIR/$METHOD/raw"
LOG_DIR="$ROOT/logs"

[[ -d "$QUERY_DIR" ]] || { echo "ERROR: query dir not found: $QUERY_DIR" >&2; exit 1; }
[[ -d "$DB_DIR"    ]] || { echo "ERROR: db dir not found: $DB_DIR" >&2; exit 1; }
mkdir -p "$RAW_DIR" "$LOG_DIR"

# ---- tools -----------------------------------------------------------------
# ssearch36 is on PATH (module fasta/36.3.8h, system default).
# phmmer / esl-sfetch come from the hmmer/3.3.2 module.
if command -v module >/dev/null 2>&1 || [[ -f /usr/share/Modules/init/bash ]]; then
  # shellcheck disable=SC1091
  source /usr/share/Modules/init/bash 2>/dev/null || true
  module load hmmer/3.3.2 2>/dev/null || true
  module load fasta 2>/dev/null || true
fi

case "$METHOD" in
  phmmer)  command -v phmmer    >/dev/null || { echo "ERROR: phmmer not found (module load hmmer/3.3.2)" >&2; exit 1; } ;;
  ssearch) command -v ssearch36 >/dev/null || { echo "ERROR: ssearch36 not found (module load fasta)" >&2; exit 1; } ;;
  *) echo "ERROR: unknown method '$METHOD' (use ssearch or phmmer)" >&2; exit 1 ;;
esac

# concatenate all queries for this group into one file (both tools iterate over
# multiple query sequences in a single invocation -> far fewer process launches)
QUERY_ALL="$RAW_DIR/.queries.fa"
cat "$QUERY_DIR"/*.fa "$QUERY_DIR"/*.faa "$QUERY_DIR"/*.fasta 2>/dev/null > "$QUERY_ALL" || true
[[ -s "$QUERY_ALL" ]] || { echo "ERROR: no query sequences found in $QUERY_DIR" >&2; exit 1; }
NQ=$(grep -c '^>' "$QUERY_ALL")

echo "==> group=$GROUP method=$METHOD evalue<$EVALUE threads=$THREADS"
echo "==> $NQ query proteins from $QUERY_DIR"

shopt -s nullglob
DBS=("$DB_DIR"/*.fasta "$DB_DIR"/*.fa "$DB_DIR"/*.faa)
[[ ${#DBS[@]} -gt 0 ]] || { echo "ERROR: no databases (*.fasta) in $DB_DIR" >&2; exit 1; }
echo "==> ${#DBS[@]} databases in $DB_DIR"

# ---- run searches ----------------------------------------------------------
for db in "${DBS[@]}"; do
  dbstem="$(basename "$db")"; dbstem="${dbstem%.*}"
  tbl="$RAW_DIR/$dbstem.tbl"
  log="$LOG_DIR/${GROUP}_${METHOD}_${dbstem}.log"

  # skip empty / non-FASTA placeholder files (databases are sometimes still
  # being staged in); a single bad db must not abort the whole run
  if [[ ! -s "$db" ]] || ! grep -q '^>' "$db"; then
    echo "    [skip] $dbstem (empty or not FASTA)"
    continue
  fi
  echo "    [$METHOD] vs $dbstem"

  # don't let one failing search kill the run under 'set -e'
  rc=0
  if [[ "$METHOD" == "phmmer" ]]; then
    # --tblout: one line per target sequence, full-sequence E-value (column 5)
    phmmer --noali --notextw --cpu "$THREADS" -E "$EVALUE" \
           --tblout "$tbl" "$QUERY_ALL" "$db" > "$log" 2>&1 || rc=$?
  else
    # -m 8: BLAST tabular (no header); E-value in column 11
    # -d 0 -b suppress full text alignments; -E sets reporting threshold
    ssearch36 -q -T "$THREADS" -E "$EVALUE" -d 0 -m 8 \
              "$QUERY_ALL" "$db" > "$tbl" 2> "$log" || rc=$?
  fi
  if [[ $rc -ne 0 ]]; then
    echo "    [warn] $METHOD failed on $dbstem (rc=$rc); see $log -- skipping"
    rm -f "$tbl"
  fi
done

# ---- collect hits into per-query multi-FASTA -------------------------------
echo "==> collecting hits (E < $EVALUE) into $OUT_DIR/*.hits.fa"
python3 "$SCRIPT_DIR/collect_hits.py" \
  --method "$METHOD" \
  --evalue "$EVALUE" \
  --db-dir "$DB_DIR" \
  --raw-dir "$RAW_DIR" \
  --query "$QUERY_ALL" \
  --out-dir "$OUT_DIR"

rm -f "$QUERY_ALL"
echo "==> done. Results in $OUT_DIR/"
