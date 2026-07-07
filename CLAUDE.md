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

By default `-g <group>` drives all three paths — queries `<group>/`, databases
`db/<group>/`, output `results/<group>/`. To mix and match (e.g. reuse one query
set against a different database collection), override individually:

```bash
-q  query group / folder        (default: -g value)
-d  database directory          (default: db/<group>; abs or root-relative path)
-o  output directory            (default: results/<group>; abs or root-relative)

# fungi queries vs the large BFD proteome set, into results/fungi_BFD/:
bin/protein_search.sh -g fungi_BFD -q fungi
```

For large database sets (e.g. `db/fungi_BFD/` ≈ 7.7k proteomes — each entry a
symlink into `shared/projects/BFD/`) do **not** run on the head node. Launch
`bin/submit_search.sh` (run it on the head node — it submits the jobs, it is not
itself `sbatch`-ed). It splits the databases into batches and submits a SLURM
job **array** (one task per batch, run in parallel), then a dependent **collect**
job that builds the `.hits.fa` files once the array succeeds:

```bash
bin/submit_search.sh -g fungi_BFD -q fungi           # 500 dbs/task, 16 concurrent
bin/submit_search.sh -b 250 -J 24 -g fungi_BFD -q fungi -m ssearch
```
`-b` sets databases per task, `-J` the max tasks running at once; all other flags
forward to `protein_search.sh`.

**Resume / caching.** Each search writes its raw table atomically and skips any
database whose table already exists, so an interrupted or partially-failed run is
resumed just by re-launching the same command — finished databases are reported
`[cached]` and only the missing ones are searched. Use `-f` to force re-running
(e.g. after changing `-e`; changing `-m` already writes to a separate `raw/` dir).

The search and collect phases can also be run separately by hand: `-s`
(search only), `-c` (collect only), `-L <dblist>` (read db paths from a file),
`-r START:COUNT` (process only that slice of the list), `-p` (print the resolved
db list and exit). These are the primitives `submit_search.sh` is built on.

### Building a tree from hits

`bin/make_tree.sh` concatenates its input FASTAs, then MUSCLE -> ClipKIT -> FastTreeMP:

```bash
# hits + an insect outgroup, explicit output prefix
bin/make_tree.sh -o results/fungi/KAN0768429.1.hits_insect \
    results/fungi/KAN0768429.1.hits.fa insect/PGA6_Apolygus_lucorum.fa

# single hits file, default prefix (first input minus its extension)
bin/make_tree.sh results/fungi/KAN0768429.1.hits.fa

# choose the amino-acid model (wag | lg | jtt; default wag)
bin/make_tree.sh -m lg results/fungi/KAN0768429.1.hits.fa
bin/make_tree.sh -h                                          # usage
```

Takes any number of input FASTAs; the tree is `<prefix>.<model>.tre` (note: this is a
FastTree ML tree, not neighbor-joining).

## Tooling / environment

UCR HPCC (Rocky 8, SLURM, environment modules). The driver auto-loads what it needs:
- `ssearch36` — module `fasta` (also on default PATH)
- `phmmer`, `esl-sfetch` — module `hmmer/3.3.2`
- `python3` with Biopython (FungiDB miniconda) — used by `bin/collect_hits.py`
- `muscle`, `clipkit`, `fasttree` — modules; used by `bin/make_tree.sh`

For long runs over large db sets, use `bin/submit_search.sh` (SLURM job array)
rather than the head node — see Commands above.

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
   followed by hits sorted best-E-value first. Hit headers are rewritten as:

   ```
   >SPECIES__HITID [db=<db_stem>] [E=<evalue>] <original description>
   ```

   `SPECIES` is taken from the `organism=` field of the source header; if absent it
   falls back to the database filename stem (with a trailing `.proteins` stripped).
   `HITID` is the original sequence id with `:;|,()` characters replaced by `_` so the
   name is safe for downstream MSA / tree tools. Output is ready for MSA.

The raw `.tbl` files are namespaced by method, but the final `<query_id>.hits.fa`
is shared per group — re-running with a different `-m` overwrites it.

3. **`bin/make_tree.sh`** — optional downstream step: concatenates one or more input
   FASTAs (e.g. a `.hits.fa` plus an outgroup), aligns with MUSCLE, trims with ClipKIT,
   and infers a maximum-likelihood tree with FastTreeMP. Generic in its inputs — see
   Commands. Outputs `<prefix>.afa`, `<prefix>.afa.clipkit`, and `<prefix>.<model>.tre`.
