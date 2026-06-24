#!/usr/bin/bash -l
#
# make_tree.sh - align a set of protein FASTAs and build a phylogenetic tree.
#
# Concatenates one or more input FASTA files, aligns with MUSCLE, trims the
# alignment with ClipKIT, and infers a tree with FastTreeMP.
#
# Usage:
#   bin/make_tree.sh [-o out_prefix] [-m model] input.fa [extra.fa ...]
#
#   -o  output prefix for all generated files
#       (default: the first input with its extension stripped)
#   -m  FastTree amino-acid model: wag | lg | jtt   (default: wag)
#   -h  show this help
#
# Examples:
#   # fungi hits + insect outgroup (the original workflow)
#   bin/make_tree.sh -o results/fungi/KAN0768429.1.hits_insect \
#       results/fungi/KAN0768429.1.hits.fa insect/PGA6_Apolygus_lucorum.fa
#
#   # just one hits file, default output prefix
#   bin/make_tree.sh results/fungi/KAN0768429.1.hits.fa
#
set -euo pipefail

OUTPREFIX=""
MODEL="wag"
while getopts "o:m:h" opt; do
  case "$opt" in
    o) OUTPREFIX="$OPTARG" ;;
    m) MODEL="$OPTARG" ;;
    h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Run '$0 -h' for usage." >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

[[ $# -ge 1 ]] || { echo "ERROR: need at least one input FASTA. '$0 -h' for usage." >&2; exit 1; }
for f in "$@"; do
  [[ -s "$f" ]] || { echo "ERROR: input not found or empty: $f" >&2; exit 1; }
done

# default output prefix derived from the first input (strip a single extension)
if [[ -z "$OUTPREFIX" ]]; then
  first="$1"; OUTPREFIX="${first%.*}"
fi
mkdir -p "$(dirname "$OUTPREFIX")"

case "$MODEL" in wag|lg|jtt) ;; *) echo "ERROR: model must be wag|lg|jtt" >&2; exit 1 ;; esac

module load fasttree
module load clipkit
module load muscle

COMBINED="$OUTPREFIX.input.fa"
ALN="$OUTPREFIX.afa"
TRIM="$ALN.clipkit"
TREE="$OUTPREFIX.${MODEL}.tre"

echo "==> combining $# input file(s) -> $COMBINED"
cat "$@" > "$COMBINED"
echo "    $(grep -c '^>' "$COMBINED") sequences"

echo "==> MUSCLE align -> $ALN"
muscle -align "$COMBINED" -output "$ALN"

echo "==> ClipKIT trim -> $TRIM"
clipkit "$ALN"        # writes $ALN.clipkit

echo "==> FastTreeMP (-$MODEL) -> $TREE"
FastTreeMP "-$MODEL" < "$TRIM" > "$TREE"

echo "==> done."
echo "    alignment: $ALN"
echo "    trimmed:   $TRIM"
echo "    tree:      $TREE"
