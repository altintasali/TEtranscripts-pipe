#!/usr/bin/env python3
"""gene_id -> gene_name lookup, extracted once from the gene GTF.

Every table this pipeline writes is keyed by `gene_id`, because that is what
TEcount/TElocal report and what all the joins depend on. With an Ensembl or
GENCODE GTF that means the outputs are full of ENSMUSG/ENSG accessions and no
readable symbol anywhere -- the symbol is sitting in the GTF's `gene_name`
attribute, unused.

This writes it out as a plain two-column table so any output can be
annotated with a join, without changing the key of a single existing file.
Deliberately a sidecar rather than an extra column on the counts matrices:
those are read with `read.delim(..., row.names = 1)` and coerced with
as.matrix (see sample_qc.R), so a character column in the middle would turn
the whole matrix character and break DESeq2 -- and on default settings
(feature_class: TE) they hold no gene rows at all, so the column would be
empty anyway.

Gzipped, like most tables this pipeline writes.  The BED tracks alongside it
in results/reference/ are the exception, not the rule -- they stay plain
because IGV cannot read plain gzip, which does not apply to a lookup table
(`zcat`/`read.delim`/`pandas.read_csv` all take .gz directly).

A GTF with no gene_name attributes (UCSC-style, where gene_id is already the
symbol, or the synthetic test fixture) still produces a complete file --
gene_name falls back to gene_id, so a join never drops rows or yields blanks
-- plus a warning, since a silently useless lookup table is worse than a
noisy one.
"""
import gzip
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_attrs(text):
    return dict(ATTR_RE.findall(text))


def main():
    gtf = snakemake.input.gtf
    out = snakemake.output[0]
    log = open(snakemake.log[0], "w")

    names = {}
    conflicts = 0
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            a = parse_attrs(cols[8])
            gid = a.get("gene_id")
            if not gid:
                continue
            name = a.get("gene_name") or a.get("gene_symbol") or ""
            prev = names.get(gid)
            if prev is None:
                names[gid] = name
            elif name and not prev:
                names[gid] = name
            elif name and prev and name != prev:
                # One gene_id carrying two different symbols is a broken
                # annotation, not something to silently pick a winner for.
                conflicts += 1

    n_named = sum(1 for v in names.values() if v)
    with gzip.open(out, "wt") as fh:
        fh.write("gene_id\tgene_name\n")
        for gid in sorted(names):
            # Fall back to the id so every row is joinable and no label is
            # ever blank downstream.
            fh.write(f"{gid}\t{names[gid] or gid}\n")

    msg = (
        f"{len(names)} genes, {n_named} with a gene_name attribute "
        f"({100.0 * n_named / len(names):.1f}%)" if names else "no genes found"
    )
    print(f"gene_id -> gene_name lookup: {msg} -> {out}", file=log)

    if names and n_named == 0:
        print(
            f"WARNING: no gene_name/gene_symbol attributes found in {gtf}. "
            "Every row falls back to gene_id, so the lookup is a no-op. This "
            "is expected for a UCSC-style GTF (gene_id is already the "
            "symbol); if you expected symbols, check that the GTF is the "
            "annotation you meant to use.",
            file=log,
        )
        print("WARNING: no gene_name attributes in the GTF; lookup is a no-op",
              file=sys.stderr)
    if conflicts:
        print(
            f"WARNING: {conflicts} feature(s) gave a gene_id a second, "
            "different gene_name; kept the first seen.",
            file=log,
        )
    log.close()


# Guarded so the module can be imported (by the unit tests) without
# running. Snakemake's script: directive executes the file with
# __name__ == "__main__", so this still runs under the workflow --
# benchmark_summary.py has been doing exactly this all along.
if __name__ == "__main__":
    main()
