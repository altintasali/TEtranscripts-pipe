#!/usr/bin/env python3
"""Classify StringTie-assembled transcripts as gene-TE chimera candidates,
using assembly structure rather than STAR chimeric-junction reads.

Complements chimera_chimeric_reads.smk's read-level chimera screen: STAR only
flags a junction as "chimeric" when a read can't be explained by one linear
(possibly spliced) alignment. A TE sitting just upstream of a gene that
splices into it via an ordinary, canonical, nearby intron aligns as a normal
spliced read and never reaches Chimeric.out.junction -- this script catches
that case instead, using StringTie's assembled transcript structure.

Method: for every multi-exon transcript in the merged StringTie GTF, exons
are ordered 5'->3' by the transcript's own strand, then, in priority order:
  1. the transcript's TSS falls inside a TE overlapping the first exon (see
     --require-tss-in-te below), and the first exon also overlaps an
     ANNOTATED transcript's first exon (--first-exons) -> a TE promoter the
     annotation already has: a chimera, but a known one, not a new
     TE-driven transcript -> annotated_promoter_embedded_te
  2. TSS-in-TE as above, a downstream exon overlaps an annotated gene exon,
     no annotated-first-exon match -> te_initiated
  3. TSS-in-TE as above, no downstream exon matches any annotated gene ->
     te_initiated_intergenic (a fully novel, TE-driven transcript)
  4. (else) last exon overlaps a TE, an earlier exon matches an annotated
     gene -> te_terminated (TE-overlapping last exon with no matching gene
     is not a real chimera -- nothing to terminate -- and is skipped);
     but if the TE (its part inside that last exon) lies inside an
     ANNOTATED transcript's last exon on the same strand (--last-exons) ->
     the gene's ordinary 3' end with a TE in its UTR ->
     annotated_terminal_exon_embedded_te (the 3' counterpart of 1.). A
     transcript extending past the annotated end into a downstream TE stays
     te_terminated.
  5. (else) an internal exon (not first, not last) overlaps a TE, and some
     other exon matches an annotated gene -> te_exonized
A transcript whose first exon overlaps a TE (and clears the TSS-in-TE gate)
never falls through to the te_terminated/te_exonized checks even if a later
exon also touches a TE -- first-exon evidence takes priority; this only
matters for the rare transcript with TEs at both ends.

--require-tss-in-te (config chimera.assembly.require_tss_in_te, default
true): by default, a first-exon TE overlap only counts toward
te_initiated/te_initiated_intergenic/annotated_promoter_embedded_te when the
transcript's own TSS (not just any part of its first exon) falls inside the
TE interval, +/- --breakpoint-tolerance. Without this, an ordinary gene with
a TE anywhere in its 5' UTR is indistinguishable from a real TE-driven
promoter. When a first-exon TE hit exists but fails the TSS gate, it is
treated exactly as if the first exon had no TE hit at all (falls through to
the te_terminated/te_exonized checks). Pass --no flag (omit the flag) to
restore the looser any-overlap behavior.

TE orientation: every classified row reports te_strand and
te_orientation_match (transcript strand vs. TE strand, "yes"/"no"/"NA") --
orientation matters for LTR-promoter chimeras, where the promoter only
drives transcription in the LTR's own orientation.

Single-exon transcripts overlapping a TE have no splice evidence to confirm
they connect to anything -> unspliced_te_only, reported separately at lower
confidence rather than classified.

Reuses the same reference tracks the junction screen already builds
(results/reference/genes.bed, exons.bed, te.bed from annotation_to_bed.py),
plus first_exons.bed (also from annotation_to_bed.py) -- no new reference
files needed.

Ambiguity: like classify_chimera_chimeric_reads.py, the reported te_id/
matched_gene_id at a multi-copy/nested locus use the first hit -- the full
overlap set is preserved in te_hits_all/gene_hits_all for inspection.

Output columns (results/chimera/assembly/transcripts.tsv):
    transcript_id, gtf_gene_id, chrom, strand, n_exons,
    transcript_start, transcript_end (the whole assembled transcript's span),
    te_exon_start, te_exon_end (the SPECIFIC exon that overlaps the TE --
        not always the first exon; use this, not transcript_start/end, to
        jump straight to the breakpoint e.g. for an IGV track),
    te_id, te_subfamily, te_family, te_class, te_strand, te_overlap_exon_rank,
    te_hits_all, te_orientation_match (transcript strand vs. TE strand),
    matched_gene_id, matched_gene_strand, gene_hits_all, strand_match,
    chimera_type

This step is ANNOTATE-ONLY like classify_chimera_chimeric_reads.py: nothing is
filtered here except -m/-c thresholds already applied by StringTie itself.
Apply your own per-sample expression / replicate-support filter downstream
(see quantify_chimera_assembly.py).

Usage:
  classify_chimera_assembly.py --gtf stringtie_merge.gtf \\
      --genes genes.bed --exons exons.bed --first-exons first_exons.bed \\
      --te te.bed \\
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
from chimera_exon_context import (
    ANNOTATED_PROMOTER,
    ANNOTATED_TERMINAL_EXON,
    ANTISENSE_TO_GENE,
)

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def load_bed(path, n_extra=0):
    """{chrom: (feats_sorted_by_start, running_max_end)} -- same interval
    index as classify_chimera_chimeric_reads.py's load_bed."""
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


def orientation_match(transcript_strand, te_strand):
    """"yes"/"no" when both strands are known, else "NA" -- used for
    te_orientation_match, the LTR-promoter-relevant signal that compares the
    assembled transcript's strand to the TE's own strand (not the matched
    gene's, which strand_match already covers)."""
    if transcript_strand not in ("+", "-") or te_strand not in ("+", "-"):
        return "NA"
    return "yes" if transcript_strand == te_strand else "no"


def order_te_hits(hits, exon_s, exon_e, strand, anchor=None, tol=0,
                  by_distance=False):
    """TE hits for one exon, the one to NAME first, by `anchor`: the
    transcript's TSS for a first exon, the last exon's splice acceptor
    (its 5' boundary) for a last exon, None for an internal exon.

    First exon: a TE containing the TSS (+/- tol), then the most bases
    inside the exon, then on the transcript's own strand, then by
    coordinate. Last exon (by_distance): the TE nearest the acceptor, then
    the same tie-breaks. Internal exon: the most bases inside the exon.

    Hits are found with +/- breakpoint tolerance, so a transcript starting
    where two TEs meet hits both -- one running through the exon and one
    lying just outside it. Taking hits[0] (lowest coordinate) named the
    outside one on "+" transcripts: on a real 84-sample run 1,250 of 2,911
    multi-TE initiation calls named a different TE than this order does,
    969 of them a TE with <= 5 bp in the first exon (e.g. an L1 ending at
    an MT2_Mm LTR promoter's TSS named instead of the MT2_Mm).

    A last exon is often a long 3' UTR holding several TEs. What makes a
    transcript TE-terminated is the gene splicing INTO a TE-derived last
    exon, so the TE nearest the splice acceptor is named -- the one the
    junction-based screens see. On the same run, renamed terminal calls
    were found by the chimeric-reads or SJ screen 41-46% of the time with
    this pick, 6-19% with the lowest coordinate, and 7-15% with the TE at
    the transcript's 3' end."""
    def at_anchor(h):
        return anchor is not None and h[0] - tol <= anchor <= h[1] + tol
    def key(h):
        inside = max(0, min(h[1], exon_e) - max(h[0], exon_s))
        if by_distance:
            dist = 0 if at_anchor(h) else min(abs(anchor - h[0]), abs(anchor - h[1]))
            return (dist, -inside, h[2][2] != strand, h[0], h[1])
        return (not at_anchor(h), -inside, h[2][2] != strand, h[0], h[1])
    return sorted(hits, key=key)


def find_gene_match(exons, exclude_rank, exons_track, chrom, tol,
                    transcript_strand="."):
    """The annotated gene this transcript's exons (transcription order,
    skipping exclude_rank) land on. Returns (gene_id, gene_strand,
    all_hit_ids) or (".", ".", []).

    A gene on the transcript's own strand wins: the first exon hitting a
    same-strand gene decides. Only when no exon hits a same-strand gene does
    the first exon with ANY gene hit decide (old behaviour), and the caller
    then types the call antisense_to_gene. Previously the first hit of any
    strand won, so a transcript through two overlapping genes on opposite
    strands could be matched to the wrong one."""
    first_any = None
    for rank, (s, e) in enumerate(exons, start=1):
        if rank == exclude_rank:
            continue
        hits = overlapping(exons_track, chrom, s - tol, e + tol)
        if not hits:
            continue
        all_ids = sorted({h[2][0] for h in hits})
        if transcript_strand in ("+", "-"):
            same = [h for h in hits if h[2][2] == transcript_strand]
            if same:
                return same[0][2][0], same[0][2][2], all_ids
        if first_any is None:
            first_any = (hits[0][2][0], hits[0][2][2], all_ids)
            if transcript_strand not in ("+", "-"):
                return first_any
    return first_any if first_any is not None else (".", ".", [])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="StringTie merged/assembled GTF")
    ap.add_argument("--genes", required=True, help="results/reference/genes.bed")
    ap.add_argument("--exons", required=True, help="results/reference/exons.bed")
    ap.add_argument("--first-exons", required=True,
                     help="results/reference/first_exons.bed")
    ap.add_argument("--last-exons", default=None,
                     help="results/reference/last_exons.bed "
                     "(annotation_splice_features.py). A TE-overlapping last "
                     "exon that overlaps an annotated last exon on the same "
                     "strand is typed annotated_terminal_exon_embedded_te "
                     "instead of te_terminated. Omitted -> not checked.")
    ap.add_argument("--te", required=True, help="results/reference/te.bed")
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument("--min-exons-for-splice-call", type=int, default=2)
    ap.add_argument("--require-tss-in-te", action="store_true",
                     help="Only call te_initiated when the transcript's own "
                          "TSS falls inside the TE, not merely its first "
                          "exon (config chimera.assembly.require_tss_in_te).")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    exons_track = load_bed(args.exons)
    first_exons_track = load_bed(args.first_exons)
    last_exons_track = load_bed(args.last_exons) if args.last_exons else None
    te = load_bed(args.te, n_extra=3)
    tol = max(args.breakpoint_tolerance, 0)

    transcripts = load_transcripts(args.gtf)

    header = [
        "transcript_id", "gtf_gene_id", "chrom", "strand", "n_exons",
        "transcript_start", "transcript_end",
        "te_exon_start", "te_exon_end",
        "te_id", "te_subfamily", "te_family", "te_class", "te_strand",
        "te_overlap_exon_rank", "te_hits_all", "te_orientation_match",
        "matched_gene_id", "matched_gene_strand", "gene_hits_all", "strand_match",
        "chimera_type",
    ]
    rows = []

    for tid, t in transcripts.items():
        ex = t["exons"]
        chrom, strand = t["chrom"], t["strand"]
        n_exons = len(ex)

        transcript_start = min(s for s, _ in ex)
        transcript_end = max(e for _, e in ex)

        if n_exons < args.min_exons_for_splice_call or strand not in ("+", "-"):
            for s, e in ex:
                hits = overlapping(te, chrom, s - tol, e + tol)
                if hits:
                    te_id, _, te_strand, te_fam, te_cls, te_sub = hits[0][2]
                    te_all = sorted({h[2][0] for h in hits})
                    rows.append([
                        tid, t["gene_id"], chrom, strand or ".", n_exons,
                        transcript_start, transcript_end, s, e,
                        te_id, te_sub, te_fam, te_cls, te_strand, 1,
                        ",".join(te_all), orientation_match(strand, te_strand),
                        ".", ".", "", "NA", "unspliced_te_only",
                    ])
                    break
            continue

        first_s, first_e = ex[0]
        last_s, last_e = ex[-1]
        te_first_hits = overlapping(te, chrom, first_s - tol, first_e + tol)
        te_last_hits = overlapping(te, chrom, last_s - tol, last_e + tol)

        # TSS = the transcript's own 5' boundary, which is first_s on "+" and
        # first_e on "-" (ex[0] is already the 5'-most exon). Gate first-exon
        # TE evidence on the TSS actually falling inside the TE, not merely
        # any overlap with the exon -- an ordinary gene with a TE anywhere in
        # its 5' UTR would otherwise also qualify. See --require-tss-in-te.
        te_first_hits_init = te_first_hits
        if args.require_tss_in_te and te_first_hits:
            tss = first_e if strand == "-" else first_s
            te_first_hits_init = [
                h for h in te_first_hits if h[0] - tol <= tss <= h[1] + tol
            ]
        te_first_hits_init = order_te_hits(
            te_first_hits_init, first_s, first_e, strand,
            anchor=first_e if strand == "-" else first_s, tol=tol)
        te_last_hits = order_te_hits(
            te_last_hits, last_s, last_e, strand,
            anchor=last_e if strand == "-" else last_s, tol=tol,
            by_distance=True)

        chimera_type = te_rank = None
        te_id = te_fam = te_cls = te_strand = "."
        te_hits_all = []
        matched_gene_id = matched_gene_strand = "."
        gene_hits_all = []
        te_exon_s = te_exon_e = None

        if te_first_hits_init:
            te_id, _, te_strand, te_fam, te_cls, te_sub = te_first_hits_init[0][2]
            te_hits_all = sorted({h[2][0] for h in te_first_hits_init})
            te_rank = 1
            te_exon_s, te_exon_e = first_s, first_e
            matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                ex, te_rank, exons_track, chrom, tol, strand
            )
            first_exon_hits = overlapping(
                first_exons_track, chrom, first_s - tol, first_e + tol
            )
            if first_exon_hits:
                # First exon lines up with a KNOWN transcript's own first
                # exon -- an annotated promoter with an embedded TE, not a
                # novel TE-driven transcript.
                chimera_type = ANNOTATED_PROMOTER
            else:
                chimera_type = "te_initiated" if matched_gene_id != "." else "te_initiated_intergenic"
        elif te_last_hits:
            te_id, _, te_strand, te_fam, te_cls, te_sub = te_last_hits[0][2]
            te_hits_all = sorted({h[2][0] for h in te_last_hits})
            te_rank = n_exons
            te_exon_s, te_exon_e = last_s, last_e
            matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                ex, te_rank, exons_track, chrom, tol, strand
            )
            if matched_gene_id != ".":
                # The 3' counterpart of annotated_promoter_embedded_te: a
                # last exon that lines up with a KNOWN transcript's own last
                # exon (same strand) is the gene's ordinary 3' end with a TE
                # in its UTR, not a new TE-terminated transcript.
                #
                # The TE itself must sit INSIDE a same-strand annotated last
                # exon (the TE portion within this transcript's last exon
                # overlaps it) -- not merely the transcript's last exon
                # overlapping one. A transcript that runs PAST the gene's
                # annotated 3' end into a downstream TE, which then supplies
                # the new end, is a real TE-terminated event and stays
                # te_terminated. (Testing exon overlap alone swept 82% of
                # te_terminated calls into this class on a real run, vs 43%
                # with the TE in an ordinary UTR.)
                annotated_last = False
                if last_exons_track is not None:
                    for ts_, te_, _ex in te_last_hits:
                        seg_s = max(ts_, last_s - tol)
                        seg_e = min(te_, last_e + tol)
                        if seg_e <= seg_s:
                            continue
                        if any(h[2][2] == strand for h in overlapping(
                                last_exons_track, chrom, seg_s, seg_e)):
                            annotated_last = True
                            break
                chimera_type = (ANNOTATED_TERMINAL_EXON if annotated_last
                                else "te_terminated")
            # else: TE-overlapping last exon but no earlier exon matches a
            # known gene -- nothing to terminate, not a real chimera; skip.
        else:
            for rank, (s, e) in enumerate(ex[1:-1], start=2):
                hits = order_te_hits(overlapping(te, chrom, s - tol, e + tol),
                                     s, e, strand)
                if hits:
                    te_id, _, te_strand, te_fam, te_cls, te_sub = hits[0][2]
                    te_hits_all = sorted({h[2][0] for h in hits})
                    te_rank = rank
                    te_exon_s, te_exon_e = s, e
                    chimera_type = "te_exonized"
                    break
            if chimera_type == "te_exonized":
                matched_gene_id, matched_gene_strand, gene_hits_all = find_gene_match(
                    ex, te_rank, exons_track, chrom, tol, strand
                )

        if chimera_type is None:
            continue

        strand_match = "NA"
        if matched_gene_strand not in (".", None):
            strand_match = "yes" if strand == matched_gene_strand else "no"

        # Strand rule: an assembled transcript on the strand opposite its
        # matched gene is antisense transcription through that gene's exon,
        # not a TE-initiated/terminated/exonized version of the gene (nor
        # the gene's own promoter). On a real 4-sample run these were the
        # te_initiated calls with the TE downstream of the whole gene and
        # the te_terminated calls with it upstream.
        if strand_match == "no" and chimera_type in (
                "te_initiated", "te_terminated", "te_exonized",
                ANNOTATED_PROMOTER, ANNOTATED_TERMINAL_EXON):
            chimera_type = ANTISENSE_TO_GENE

        rows.append([
            tid, t["gene_id"], chrom, strand, n_exons,
            transcript_start, transcript_end, te_exon_s, te_exon_e,
            te_id, te_sub, te_fam, te_cls, te_strand, te_rank,
            ",".join(te_hits_all), orientation_match(strand, te_strand),
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
