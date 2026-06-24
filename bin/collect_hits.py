#!/usr/bin/env python3
"""Collect significant search hits into one multi-FASTA per query protein.

Reads the raw tabular output produced by protein_search.sh (phmmer --tblout or
SSEARCH -m 8), filters by E-value, fetches each hit sequence from its source
database, and writes results/<group>/<query_id>.hits.fa.

Each output record header is annotated with its source database and E-value:
    >TARGET_ID [db=<database_stem>] [E=<evalue>] <original description>

The query protein itself is written as the first record so the file is ready
for downstream multiple-sequence alignment.
"""
import argparse
import glob
import os
import re
import sys

from Bio import SeqIO


def parse_tbl(path, method, threshold):
    """Yield (query_id, target_id, evalue) for hits with evalue < threshold."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            try:
                if method == "phmmer":
                    # tblout: target(0) ... query(2) ... full-seq E-value(4)
                    target, query, ev = f[0], f[2], float(f[4])
                else:  # ssearch -m 8: query(0) target(1) ... E-value(10)
                    query, target, ev = f[0], f[1], float(f[10])
            except (IndexError, ValueError):
                continue
            if ev < threshold:
                yield query, target, ev


def safe_name(name):
    """Filesystem-safe query id for the output filename."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=["phmmer", "ssearch"])
    ap.add_argument("--evalue", required=True, type=float)
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--query", required=True, help="concatenated query FASTA")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    # index query sequences (for the leading record in each output file)
    queries = SeqIO.index(args.query, "fasta")

    # lazily index each database FASTA for fast random access by id
    db_index = {}

    def get_db(stem):
        if stem not in db_index:
            for ext in (".fasta", ".fa", ".faa"):
                p = os.path.join(args.db_dir, stem + ext)
                if os.path.exists(p):
                    db_index[stem] = SeqIO.index(p, "fasta")
                    break
            else:
                db_index[stem] = None
        return db_index[stem]

    # gather hits per query: {query_id: [(evalue, db_stem, target_id)]}
    hits = {}
    n_tables = 0
    for tbl in sorted(glob.glob(os.path.join(args.raw_dir, "*.tbl"))):
        n_tables += 1
        db_stem = os.path.splitext(os.path.basename(tbl))[0]
        for query, target, ev in parse_tbl(tbl, args.method, args.evalue):
            hits.setdefault(query, []).append((ev, db_stem, target))

    if n_tables == 0:
        sys.exit("ERROR: no raw .tbl files found in " + args.raw_dir)

    os.makedirs(args.out_dir, exist_ok=True)
    total_hits = 0
    for query in sorted(hits):
        out = os.path.join(args.out_dir, safe_name(query) + ".hits.fa")
        seen = set()
        with open(out, "w") as oh:
            # query first, clearly labelled
            if query in queries:
                qrec = queries[query]
                oh.write(f">{qrec.id} [QUERY] {qrec.description[len(qrec.id):].strip()}\n")
                oh.write(str(qrec.seq) + "\n")
            # hits, best E-value first
            for ev, db_stem, target in sorted(hits[query], key=lambda x: x[0]):
                key = (db_stem, target)
                if key in seen:
                    continue
                seen.add(key)
                db = get_db(db_stem)
                if db is None or target not in db:
                    sys.stderr.write(f"  warn: {target} not found in db {db_stem}\n")
                    continue
                rec = db[target]
                desc = rec.description[len(rec.id):].strip()
                # prefix the hit id with the species name from the organism= field
                m = re.search(r"organism=(\S+)", desc)
                species = m.group(1) if m else db_stem
                # sanitize ':' (and ';' '|' ',' '(' ')') in the id so the name is
                # safe for downstream MSA / tree-building tools
                hit_id = re.sub(r"[:;|,()]", "_", rec.id)
                oh.write(f">{species}__{hit_id} [db={db_stem}] [E={ev:.1e}] {desc}\n")
                oh.write(str(rec.seq) + "\n")
                total_hits += 1
        print(f"    {os.path.basename(out)}: {len(seen)} hits")

    print(f"==> {total_hits} hits across {len(hits)} query proteins")


if __name__ == "__main__":
    main()
