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
#   bin/protein_search.sh [-g group] [-q query_group] [-d db_dir] [-o out_dir] \
#                         [-m method] [-e evalue] [-t threads] \
#                         [-L dblist] [-r START:COUNT] [-s] [-c] [-f] [-p]
#
#   -g  major group name; sets the defaults for -q, -d and -o   (default: fungi)
#   -q  query group / folder name                               (default: -g value)
#   -d  database directory                                      (default: db/<group>)
#   -o  output directory                                        (default: results/<group>)
#   -m  search method: ssearch | phmmer                         (default: phmmer)
#   -e  E-value threshold (hits must be <)                      (default: 1e-20)
#   -t  CPU threads per search                                  (default: 8)
#   -L  read database list from this file (one path per line) instead of scanning db_dir
#   -r  process only COUNT databases starting at line START (1-based) of the list
#   -s  search only (write raw tables, skip the collect step)
#   -c  collect only (skip searches, just (re)build the .hits.fa files)
#   -f  force: re-run searches even if a cached raw table already exists
#   -p  print the resolved database list to stdout and exit (no search/collect)
#
# -d and -o accept absolute paths or paths relative to the project root.
#
# Resume / caching: a search writes its raw table atomically and, by default,
# skips any database whose table already exists. So an interrupted run is resumed
# simply by launching it again. Use -f to ignore the cache (e.g. after changing
# -e); changing -m already writes to a separate raw/ directory.
#
# Examples:
#   bin/protein_search.sh                          # fungi, phmmer, 1e-20
#   bin/protein_search.sh -g fungi -m ssearch      # Smith-Waterman via SSEARCH
#   bin/protein_search.sh -g insect -e 1e-30 -t 16
#   # fungi queries vs a different db set, into results/fungi_BFD:
#   bin/protein_search.sh -g fungi_BFD -q fungi
#   # search one 500-db slice only, then collect separately (see submit_search.sh):
#   bin/protein_search.sh -g fungi_BFD -q fungi -L dblist.txt -r 1:500 -s
#   bin/protein_search.sh -g fungi_BFD -q fungi -c
#
set -euo pipefail

# ---- defaults --------------------------------------------------------------
GROUP="fungi"
QUERY_GROUP=""   # default filled from GROUP below
DB_DIR=""        # default filled from GROUP below
OUT_DIR=""       # default filled from GROUP below
METHOD="phmmer"
EVALUE="1e-20"
THREADS="8"
DBLIST=""        # -L: explicit manifest of db paths
RANGE=""         # -r: START:COUNT slice of the db list
DO_SEARCH=1
DO_COLLECT=1
FORCE=0
PRINT_DBS=0

while getopts "g:q:d:o:m:e:t:L:r:scfph" opt; do
  case "$opt" in
    g) GROUP="$OPTARG" ;;
    q) QUERY_GROUP="$OPTARG" ;;
    d) DB_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    m) METHOD="$OPTARG" ;;
    e) EVALUE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    L) DBLIST="$OPTARG" ;;
    r) RANGE="$OPTARG" ;;
    s) DO_COLLECT=0 ;;
    c) DO_SEARCH=0 ;;
    f) FORCE=1 ;;
    p) PRINT_DBS=1 ;;
    h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Run '$0 -h' for usage." >&2; exit 1 ;;
  esac
done

# ---- resolve project root (parent of bin/) --------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

# fill in defaults from GROUP for any path not overridden on the command line
: "${QUERY_GROUP:=$GROUP}"
: "${DB_DIR:=db/$GROUP}"
: "${OUT_DIR:=results/$GROUP}"

# allow -d / -o to be relative to the project root or absolute
abspath() { case "$1" in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$ROOT" "$1" ;; esac; }
QUERY_DIR="$ROOT/$QUERY_GROUP"
DB_DIR="$(abspath "$DB_DIR")"
OUT_DIR="$(abspath "$OUT_DIR")"
RAW_DIR="$OUT_DIR/$METHOD/raw"
LOG_DIR="$ROOT/logs"

[[ -d "$DB_DIR" ]] || { echo "ERROR: db dir not found: $DB_DIR" >&2; exit 1; }

# ---- enumerate databases ---------------------------------------------------
# (db entries may be symlinks, so do NOT filter with -type f; broken/empty links
#  are skipped at search time by the per-db validity check below)
if [[ -n "$DBLIST" ]]; then
  mapfile -t DBS < "$(abspath "$DBLIST")"
else
  mapfile -t DBS < <(find "$DB_DIR" -maxdepth 1 \
      \( -name '*.fasta' -o -name '*.fa' -o -name '*.faa' \) | sort)
fi

# -r START:COUNT -> keep COUNT entries starting at 1-based line START
if [[ -n "$RANGE" ]]; then
  start="${RANGE%%:*}"; count="${RANGE##*:}"
  DBS=("${DBS[@]:start-1:count}")
fi

[[ ${#DBS[@]} -gt 0 ]] || { echo "ERROR: no databases found (db_dir=$DB_DIR list=${DBLIST:-none} range=${RANGE:-all})" >&2; exit 1; }

if [[ $PRINT_DBS -eq 1 ]]; then
  printf '%s\n' "${DBS[@]}"
  exit 0
fi

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

if [[ $DO_SEARCH -eq 1 ]]; then
  case "$METHOD" in
    phmmer)  command -v phmmer    >/dev/null || { echo "ERROR: phmmer not found (module load hmmer/3.3.2)" >&2; exit 1; } ;;
    ssearch) command -v ssearch36 >/dev/null || { echo "ERROR: ssearch36 not found (module load fasta)" >&2; exit 1; } ;;
    *) echo "ERROR: unknown method '$METHOD' (use ssearch or phmmer)" >&2; exit 1 ;;
  esac
fi

# ---- queries (needed by both search and collect) ---------------------------
# unique temp so parallel array tasks never clobber a shared file
[[ -d "$QUERY_DIR" ]] || { echo "ERROR: query dir not found: $QUERY_DIR" >&2; exit 1; }
QUERY_ALL="$(mktemp "${TMPDIR:-/tmp}/queries.XXXXXX.fa")"
trap 'rm -f "$QUERY_ALL"' EXIT
cat "$QUERY_DIR"/*.fa "$QUERY_DIR"/*.faa "$QUERY_DIR"/*.fasta 2>/dev/null > "$QUERY_ALL" || true
[[ -s "$QUERY_ALL" ]] || { echo "ERROR: no query sequences found in $QUERY_DIR" >&2; exit 1; }
NQ=$(grep -c '^>' "$QUERY_ALL")

echo "==> group=$GROUP query=$QUERY_GROUP method=$METHOD evalue<$EVALUE threads=$THREADS"
echo "==> $NQ query proteins from $QUERY_DIR; ${#DBS[@]} databases to process"

# ---- run searches ----------------------------------------------------------
if [[ $DO_SEARCH -eq 1 ]]; then
  for db in "${DBS[@]}"; do
    dbstem="$(basename "$db")"; dbstem="${dbstem%.*}"
    tbl="$RAW_DIR/$dbstem.tbl"
    log="$LOG_DIR/${GROUP}_${METHOD}_${dbstem}.log"

    # resume: a finished search renamed its table into place, so its presence
    # means "already done" (even a 0-hit search produces a valid table)
    if [[ $FORCE -eq 0 && -e "$tbl" ]]; then
      echo "    [cached] $dbstem"
      continue
    fi

    # skip empty / broken-symlink / non-FASTA entries; one bad db must not abort
    if [[ ! -s "$db" ]] || ! grep -q '^>' "$db"; then
      echo "    [skip] $dbstem (empty, missing, or not FASTA)"
      continue
    fi
    echo "    [$METHOD] vs $dbstem"

    # write to a temp table, rename only on success -> never leave a partial
    # table that a resumed run would mistake for a completed search
    tmp="$tbl.partial.$$"
    rc=0
    if [[ "$METHOD" == "phmmer" ]]; then
      # --tblout: one line per target sequence, full-sequence E-value (column 5)
      phmmer --noali --notextw --cpu "$THREADS" -E "$EVALUE" \
             --tblout "$tmp" "$QUERY_ALL" "$db" > "$log" 2>&1 || rc=$?
    else
      # -m 8: BLAST tabular (no header); E-value in column 11
      ssearch36 -q -T "$THREADS" -E "$EVALUE" -d 0 -m 8 \
                "$QUERY_ALL" "$db" > "$tmp" 2> "$log" || rc=$?
    fi
    if [[ $rc -eq 0 ]]; then
      mv -f "$tmp" "$tbl"
    else
      echo "    [warn] $METHOD failed on $dbstem (rc=$rc); see $log -- skipping"
      rm -f "$tmp"
    fi
  done
fi

# ---- collect hits into per-query multi-FASTA -------------------------------
if [[ $DO_COLLECT -eq 1 ]]; then
  echo "==> collecting hits (E < $EVALUE) into $OUT_DIR/*.hits.fa"
  python3 "$SCRIPT_DIR/collect_hits.py" \
    --method "$METHOD" \
    --evalue "$EVALUE" \
    --db-dir "$DB_DIR" \
    --raw-dir "$RAW_DIR" \
    --query "$QUERY_ALL" \
    --out-dir "$OUT_DIR"
  echo "==> done. Results in $OUT_DIR/"
else
  echo "==> search phase done ($RAW_DIR). Run with -c to collect."
fi
