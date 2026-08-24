#!/usr/bin/env python3
"""Classify StringTie-assembled transcripts as gene-TE chimera candidates,
using assembly structure rather than STAR chimeric-junction reads.

Complements chimera_junction.smk's read-level chimera screen: STAR only
flags a junction as "chimeric" when a read can't be explained by one linear
(possibly spliced) alignment. A TE sitting just upstream of a gene that
splices into it via an ordinary, canonical, nearby intron aligns as a normal
spliced read and never reaches Chimeric.out.junction -- this script catches
that case instead, using StringTie's assembled transcript structure.

Method: for every multi-exon transcript in the merged StringTie GTF, exons
are ordered 5'->3' by the transcript's own strand, then, in priority order:
  1. first exon overlaps a TE, a downstream exon overlaps an annotated gene
     exon -> te_initiated
  2. first exon overlaps a TE, no downstream exon matches any annotated gene
     -> te_initiated_intergenic (a fully novel, TE-driven transcript)
  3. (else) last exon overlaps a TE, an earlier exon matches an annotated
     gene -> te_terminated (TE-overlapping last exon with no matching gene
     is not a real chimera -- nothing to terminate -- and is skipped)
  4. (else) an internal exon (not first, not last) overlaps a TE, and some
     other exon matches an annotated gene -> te_exonized
A transcript whose first exon overlaps a TE never falls through to the
te_terminated/te_exonized checks even if a later exon also touches a TE --
first-exon evidence takes priority; this only matters for the rare
transcript with TEs at both ends.

Single-exon transcripts overlapping a TE have no splice evidence to confirm
they connect to anything -> unspliced_te_only, reported separately at lower
confidence rather than classified.

Reuses the same reference tracks the junction screen already builds
(results/reference/genes.bed, exons.bed, te.bed from annotation_to_bed.py) --
no new reference files needed.

Ambiguity: like parse_chimeric_junctions.py, the reported te_id/
matched_gene_id at a multi-copy/nested locus use the first hit -- the full
overlap set is preserved in te_hits_all/gene_hits_all for inspection.

Output columns (results/chimera_assembly/candidates.tsv):
    transcript_id, gtf_gene_id, chrom, strand, n_exons,
    first_exon_start, first_exon_end,
    te_id, te_family, te_class, te_overlap_exon_rank, te_hits_all,
    matched_gene_id, matched_gene_strand, gene_hits_all, strand_match,
    chimera_type

This step is ANNOTATE-ONLY like parse_chimeric_junctions.py: nothing is
filtered here except -m/-c thresholds already applied by StringTie itself.
Apply your own per-sample expression / replicate-support filter downstream
(see quantify_chimera_assembly.py).

Usage:
  classify_chimera_assembly.py --gtf stringtie_merge.gtf \\
      --genes genes.bed --exons exons.bed --te te.bed \\
      --breakpoint-tolerance 5 --out candidates.tsv.gz
"""
import argparse
import bisect
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_write

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def load_bed(path, n_extra=0):
    """{chrom: (feats_sorted_by_start, running_max_end)} -- same interval
    index as parse_chimeric_junctions.py's load_bed."""
    raw = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            extras = tuple(cols[3 : 6 + n_extra])
            raw.setdefault(chrom, []).append((start, end, extras))
    tracks = {}
    for chrom, feats in raw.items():
        feats.sort()
        running = float("-inf")
        max_end = []
        for _, e, _ in feats:
            running = max(running, e)
            max_end.append(running)
        tracks[chrom] = (feats, max_end)
    return tracks


def overlapping(track, chrom, start0, end0):
    if chrom not in track:
        return []
    feats, max_end = track[chrom]
    lo = bisect.bisect_right(max_end, start0)
    hits = []
    for i in range(lo, len(feats)):
        s, e, extras = feats[i]
        if s >= end0:
            break
        if e > start0:
            hits.append((s, e, extras))
    return hits


def load_transcripts(gtf_path):
    """{transcript_id: {"chrom", "strand", "gene_id", "exons": [(s, e), ...]}}
    Exons are returned in transcription order (5' -> 3'), using strand."""
    transcripts = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom, start, end, strand = cols[0], int(cols[3]) - 1, int(cols[4]), cols[6]
            attrs = parse_attrs(cols[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            t = transcripts.setdefault(
                tid, {"chrom": chrom, "strand": strand,
                      "gene_id": attrs.get("gene_id", "."), "exons": []}
            )
            t["exons"].append((start, end))
    for t in transcripts.values():
        t["exons"].sort()
        if t["strand"] == "-":
            t["exons"].reverse()
    return transcripts


def find_gene_match(exons, exclude_rank, exons_track, chrom, tol):
    """First exon (transcription-order rank, skipping exclude_rank) that
    overlaps an annotated gene exon. Returns (gene_id, gene_strand,
    all_hit_ids) or (".", ".", [])."""
    for rank, (s, e) in enumerate(exons, start=1):
        if rank == exclude_rank:
            continue
        hits = overlapping(exons_track, chrom, s - tol, e + tol)
        if hits:
            all_ids = sorted({h[2][0] for h in hits})
            return hits[0][2][0], hits[0][2][2], all_ids
    return ".", ".", []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="StringTie merged/assembled GTF")
    ap.add_argument("--genes", required=True, help="results/reference/genes.bed")
    ap.add_argument("--exons", required=True, help="results/reference/exons.bed")
    ap.add_argument("--te", required=True, help="results/reference/te.bed")
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument("--min-exons-for-splice-call", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    exons_track = load_bed(args.exons)
    te = load_bed(args.te, n_extra=2)
    tol = max(args.breakpoint_tolerance, 0)

    transcripts = load_transcripts(args.gtf)

    header = [
        "transcript_id", "gtf_gene_id", "chrom", "strand", "n_exons",
        "first_exon_start", "first_exon_end",
        "te_id", "te_family", "te_class", "te_overlap_exon_rank", "te_hits_all",
        "matched_gene_id", "matched_gene_strand", "gene_hits_all", "strand_match",
        "chimera_type",
    ]
    rows = []

    for tid, t in transcripts.items():
        ex = t["exons"]
        chrom, strand = t["chrom"], t["strand"]
        n_exons = len(ex)

        if n_exons < args.min_exons_for_splice_call or strand not in ("+", "-"):
            for s, e in ex:
                hits = overlapping(te, chrom, s - tol, e + tol)
                if hits:
                    te_id, _, te_strand, te_fam, te_cls = hits[0][2]
                    te_all = sorted({h[2][0] for h in hits})
                    rows.append([
                        tid, t["gene_id"], chrom, strand or ".", n_exons,
                        s, e, te_id, te_fam, te_cls, 1, ",".join(te_all),
                        ".", ".", "", "NA", "unspliced_te_only",
                    ])
                    break
            continue

        first_s, first_e = ex[0]
        last_s, last_e = ex[-1]
        te_first_hits = overlapping(te, chrom, first_s - tol, first_e + tol)
        te_last_hits = overlapping(te, chrom, last_s - tol, last_e + tol)

        chimera_type = te_rank = None
        te_id = te_fam = te_cls = "."
        te_hits_all = []
        matched_gene_id = matched_gene_strand = "."
        gene_hits_all = []

        if te_first_hits:
            te_id, _, te_strand, te_fam, te_cls = te_first_hits[0][2]
            te_hits_all = sorted({h[2][0] for h in te_first_hits})
            te_rank = 1
            matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                ex, te_rank, exons_track, chrom, tol
            )
            chimera_type = "te_initiated" if matched_gene_id != "." else "te_initiated_intergenic"
        elif te_last_hits:
            te_id, _, te_strand, te_fam, te_cls = te_last_hits[0][2]
            te_hits_all = sorted({h[2][0] for h in te_last_hits})
            te_rank = n_exons
            matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                ex, te_rank, exons_track, chrom, tol
            )
            if matched_gene_id != ".":
                chimera_type = "te_terminated"
            # else: TE-overlapping last exon but no earlier exon matches a
            # known gene -- nothing to terminate, not a real chimera; skip.
        else:
            for rank, (s, e) in enumerate(ex[1:-1], start=2):
                hits = overlapping(te, chrom, s - tol, e + tol)
                if hits:
                    te_id, _, te_strand, te_fam, te_cls = hits[0][2]
                    te_hits_all = sorted({h[2][0] for h in hits})
                    te_rank = rank
                    chimera_type = "te_exonized"
                    break
            if chimera_type == "te_exonized":
                matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                    ex, te_rank, exons_track, chrom, tol
                )

        if chimera_type is None:
            continue

        strand_match = "NA"
        if matched_gene_strand not in (".", None):
            strand_match = "yes" if strand == matched_gene_strand else "no"

        rows.append([
            tid, t["gene_id"], chrom, strand, n_exons,
            first_s, first_e, te_id, te_fam, te_cls, te_rank, ",".join(te_hits_all),
            matched_gene_id, matched_gene_strand, ",".join(gene_hits_all), strand_match,
            chimera_type,
        ])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{len(rows)} TE-related transcript events from {args.gtf}")


if __name__ == "__main__":
    main()
