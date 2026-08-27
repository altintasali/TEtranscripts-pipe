#!/usr/bin/env python3
"""Convert the reference GTF files into the BED tracks the chimera junction
screen overlaps breakpoints against.

Three outputs (all BED, sorted, in --outdir):

  genes.bed   BED6 per gene:      chrom, start, end, gene_id, score, strand
  exons.bed   BED6 per exon:      chrom, start, end, gene_id, score, strand
  te.bed      BED9 per TE INSERTION:
                                  chrom, start, end, te_id, score, strand,
                                  family, class, subfamily
              te_id is the individual copy (transcript_id, e.g. L1PA2_dup1);
              subfamily is gene_id (e.g. L1PA2) -- the level TEcount reports,
              kept so chimera calls can still be tied back to a TEcount row.
              Note subfamily != family (L1PA2 vs L1).

Gene exons come from the annotation GTF (ref.gtf, the same file STAR's index
uses). TE insertions come from ref.te_gtf -- the curated TE GTF (TEtranscripts
rmsk-based format, `exon` features carrying gene_id/family_id/class_id
attributes), i.e. the same input the user already maintains. Nothing here
depends on the library being stranded or paired-end: the tracks carry strand
annotations, but the breakpoint-overlap test itself is strand-agnostic (strand
is recorded per junction by the parser for downstream use).

Usage:
    annotation_to_bed.py --gtf GENE.GTF --te-gtf TE.GTF --outdir OUTDIR
"""
import argparse
import os
import sys


def parse_attrs(attr_text):
    """GTF attribute field -> {key: value} (values unquoted)."""
    attrs = {}
    for chunk in attr_text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(None, 1)
        if len(parts) == 2:
            key, value = parts
            attrs[key] = value.strip('"')
        elif len(parts) == 1:
            attrs[parts[0]] = ""
    return attrs


def emit(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--te-gtf", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if not os.path.isfile(args.gtf) or not os.path.isfile(args.te_gtf):
        sys.exit(f"error: --gtf / --te-gtf must point at existing files "
                 f"(got {args.gtf!r}, {args.te_gtf!r})")
    os.makedirs(args.outdir, exist_ok=True)

    # gene_id -> {strand, chrom, start, end, family, class}
    genes = {}
    exons = []
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = (
                cols[0], cols[1], cols[2], cols[3], cols[4], cols[5],
                cols[6], cols[7], cols[8],
            )
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            gid = a.get("gene_id")
            if not gid:
                continue
            g = genes.setdefault(
                gid, {"chrom": chrom, "strand": strand, "start": s, "end": e}
            )
            if chrom != g["chrom"]:
                g["chrom"] = chrom  # keep first; multi-contig genes are pathological
            g["start"] = min(g["start"], s)
            g["end"] = max(g["end"], e)
            if feat == "exon":
                exons.append((chrom, s, e, gid, strand))

    # TE GTF: features are the insertion loci themselves (exon feature in the
    # TEtranscripts rmsk-derived GTFs); keep subfamily/family/class for the
    # report.
    #
    # Keyed by (chrom, transcript_id) -- the individual insertion -- NOT by
    # gene_id.  In the TEtranscripts curated GTF convention gene_id is the
    # SUBFAMILY (L1PA2) and transcript_id the individual copy (L1PA2_dup1),
    # the same convention build_telocal_index.py relies on.  Keying by
    # gene_id merged every copy of a subfamily into one min/max bounding
    # box: on a real annotation that produced single intervals spanning
    # ~192 Mb (whole chromosomes), so the breakpoint-overlap test could not
    # fail and essentially every junction was called a TE hit.  Multiple
    # rows sharing one transcript_id are still merged -- that is a
    # fragmented single insertion, where min/max IS correct.
    te_loci = {}
    with open(args.te_gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = cols[:9]
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            # fall back to gene_id for non-TEtranscripts TE GTFs that carry
            # no transcript_id (then it degrades to the old behaviour rather
            # than dropping the feature entirely)
            tid = a.get("transcript_id") or a.get("gene_id")
            if not tid:
                continue
            # A real TE GTF has millions of insertions, so this dict is the
            # rule's whole memory footprint. Store a plain list (far lighter
            # than a per-row dict) and intern the handful of repeated
            # strings -- chrom, strand, subfamily, family, class have a few
            # dozen to ~1k distinct values across millions of rows, so
            # interning collapses them to one object each. Only tid is
            # genuinely unique per row.
            key = (sys.intern(chrom), tid)
            g = te_loci.get(key)
            if g is None:
                te_loci[key] = [
                    s, e, sys.intern(strand),
                    sys.intern(a.get("gene_id", ".")),
                    sys.intern(a.get("family_id", ".")),
                    sys.intern(a.get("class_id", ".")),
                ]
            else:
                if s < g[0]:
                    g[0] = s
                if e > g[1]:
                    g[1] = e

    gene_rows = sorted(
        (g["chrom"], g["start"] - 1, g["end"], gid, ".", g["strand"])
        for gid, g in genes.items()
    )
    exon_rows = sorted(
        (chrom, s - 1, e, gid, ".", strand) for chrom, s, e, gid, strand in exons
    )
    # Sort the keys in place and write directly: building a second list of
    # 3.7M formatted tuples doubled peak RSS for no benefit.
    # key mirrors the old (chrom, start, end, name) ordering exactly, so
    # te.bed stays byte-identical to what the previous implementation
    # produced rather than merely equivalent.
    te_order = sorted(
        te_loci, key=lambda k: (k[0], te_loci[k][0], te_loci[k][1], k[1])
    )
    n_te = len(te_order)

    emit(os.path.join(args.outdir, "genes.bed"), gene_rows)
    emit(os.path.join(args.outdir, "exons.bed"), exon_rows)
    with open(os.path.join(args.outdir, "te.bed"), "w") as fh:
        for k in te_order:
            chrom, tid = k
            start, end, strand, subfamily, family, cls = te_loci[k]
            fh.write(
                f"{chrom}\t{start - 1}\t{end}\t{tid}\t.\t{strand}\t"
                f"{family}\t{cls}\t{subfamily}\n"
            )
    print(
        f"wrote {len(gene_rows)} genes, {len(exon_rows)} exons, "
        f"{n_te} TE loci to {args.outdir}"
    )


if __name__ == "__main__":
    main()
