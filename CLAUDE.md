# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Homology search project: take query proteins (currently polygalacturonases, "PGA")
and find significant homologs across grouped reference proteome databases, producing
one alignment-ready multi-FASTA per query protein.

## Layout

```
<group>/                 query proteins for a major taxonomic group
  fungi/PGA_Ndisc.fa       (e.g. fungi/, insect/) — FASTA, one or more proteins
db/<group>/              subject databases for that group
  db/fungi/*.fasta         FungiDB-68 annotated proteomes, one file per genome
results/<group>/         OUTPUT: <query_id>.hits.fa per query protein
  results/<group>/<method>/raw/*.tbl   raw per-database search tables
bin/                     the search framework
logs/                    per-database search logs
```

The group name is shared between query folder and db folder: `fungi/` queries are
searched against every database in `db/fungi/`. Add a new group by creating `<group>/`
(queries) and `db/<group>/` (databases) — no code changes needed.

## Commands

```bash
# Default: fungi group, phmmer, E < 1e-20, 8 threads
bin/protein_search.sh

# Choose group / method / threshold / threads
bin/protein_search.sh -g fungi  -m phmmer  -e 1e-20 -t 8
bin/protein_search.sh -g fungi  -m ssearch -e 1e-20 -t 8     # Smith-Waterman
bin/protein_search.sh -g insect -e 1e-30  -t 16
bin/protein_search.sh -h                                     # usage
```

`-m phmmer` (default) uses HMMER profile-vs-sequence search; `-m ssearch` uses
full Smith-Waterman (FASTA `ssearch36`). They differ slightly in sensitivity.

## Tooling / environment

UCR HPCC (Rocky 8, SLURM, environment modules). The driver auto-loads what it needs:
- `ssearch36` — module `fasta` (also on default PATH)
- `phmmer`, `esl-sfetch` — module `hmmer/3.3.2`
- `python3` with Biopython (FungiDB miniconda) — used by `bin/collect_hits.py`

For long runs over large db sets, submit via `sbatch` rather than the head node.

## Architecture

Two-stage pipeline:

1. **`bin/protein_search.sh`** — concatenates all queries in `<group>/` into one
   file (both tools iterate over multiple queries per invocation, minimizing process
   launches), then runs the chosen method once per database, writing a raw table to
   `results/<group>/<method>/raw/<db_stem>.tbl`. The `-E` threshold is passed to the
   search tool so only significant hits are reported.

2. **`bin/collect_hits.py`** — parses every raw table (phmmer `--tblout` full-sequence
   E-value, column 5; or ssearch `-m 8` E-value, column 11), filters `evalue < threshold`,
   fetches each hit sequence from its source database via `Bio.SeqIO.index`, and writes
   `results/<group>/<query_id>.hits.fa`. Each file leads with the query (`[QUERY]`)
   followed by hits sorted best-E-value first; hit headers are annotated
   `[db=<genome>] [E=<evalue>]`. Output is ready for MSA.

The raw `.tbl` files are namespaced by method, but the final `<query_id>.hits.fa`
is shared per group — re-running with a different `-m` overwrites it.
