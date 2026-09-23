#!/usr/bin/env python3
"""Classify STAR's normal splice junctions (SJ.out.tab) as gene-TE chimera
candidates -- a third, assembly-free line of evidence alongside
classify_chimera_reads.py (STAR chimeric-junction reads) and
classify_chimera_assembly.py (StringTie assembly structure).

Why this exists: chimera_reads.smk's screen only sees a junction when STAR
cannot explain the read as one linear (possibly spliced) alignment.  A TE
sitting just upstream of a gene, spliced into it through a normal, canonical,
nearby intron, aligns as a completely ordinary spliced read and never reaches
Chimeric.out.junction at all -- chimera_assembly.smk's StringTie screen exists
to catch exactly that same blind spot, but needs enough per-sample coverage to
successfully ASSEMBLE a multi-exon transcript.  This screen catches the same
blind spot at the individual READ-JUNCTION level instead: SJ.out.tab already
carries one row per distinct splice junction found anywhere in the sample,
independent of whether StringTie could assemble anything from it.

Unlike Chimeric.out.junction (one row per chimeric READ, many reads sharing
one junction, requiring per-sample accumulation), STAR's SJ.out.tab already
reports one deduplicated row per distinct (chrom, intron_start, intron_end,
strand) junction, with its own unique_reads/multi_reads counts -- so this
script needs no read-level grouping, unlike classify_chimera_reads.py.

SJ.out.tab is always CIS (same-chromosome) by construction -- a normal
spliced alignment cannot jump chromosomes -- so, unlike the chimeric-read
screen, chimera_type here is never "." for a gene<->TE event with a known
strand.

Column format (STAR manual, and merge_splice_junctions.py's reference
parsing): chrom, intron_start (1-based, first intron base), intron_end
(1-based, last intron base), strand_code (0=undefined, 1=+, 2=-), motif
(0=non-canonical, 1-6=canonical GT/AG family), annotated (0=novel,
1=already in the annotation), unique_reads, multi_reads, max_overhang.

Donor/acceptor mapping: the exon ending at intron_start-1 and the exon
starting at intron_end+1 are the two splice-site breakpoints. Which one is
the transcript's 5' "donor" site vs. 3' "acceptor" site depends on strand
(on "-", transcription runs from the higher genomic coordinate to the
lower one, so intron_end+1 is the donor side and intron_start-1 the
acceptor side). When strand is undefined ("."), the "+"-strand convention
is used as an arbitrary but consistent default -- in practice this only
affects non-canonical (motif 0) junctions, which --require-canonical
excludes by default anyway (see below).

--require-canonical (default ON here, opposite of chimera_reads's default
require_canonical_junction: false): unlike a genuinely chimeric read, this
screen's whole value proposition -- working on unstranded data via the
splice motif -- depends on STAR having assigned a real strand, which only
happens for canonical (motif != 0) junctions. Non-canonical "junctions" in
SJ.out.tab are also disproportionately alignment noise on their own.

--min-unique-reads: STAR's own per-junction quality signal (uniquely
mapping reads spanning the junction) -- a lightweight floor against
single-read noise, independent of --require-canonical.

Output columns (results/chimera/sj/per_sample/{sample}_te_gene_junctions.tsv):
    event_id, sample, chrom, intron_start, intron_end, strand, motif,
    canonical, annotated, unique_reads, multi_reads, overhang,
    donor_hits, acceptor_hits, direction, direction_ambiguous,
    gene_id, gene_strand, te_id, te_subfamily, te_family, te_class,
    chimera_type, antisense_flag, library_strand, transcript_strand,
    gene_strand_match

This step is ANNOTATE-ONLY like its two siblings: --min-unique-reads and
--require-canonical are STAR-native quality gates, not expression/evidence
filters. Apply your own filter downstream.

Usage:
  classify_chimera_sj.py --sj SJ.out.tab --genes genes.bed --exons exons.bed \\
      --te te.bed --sample S1 --out S1_te_gene_junctions.tsv.gz
"""
import argparse
import bisect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_write


def load_bed(path, n_extra=0):
    """{chrom: (feats_sorted_by_start, running_max_end)} -- same interval
    index as classify_chimera_reads.py's load_bed."""
    raw = {}
    with open(path) as fh:
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


CANONICAL_MOTIFS = {1, 2, 3, 4, 5, 6}  # anything but 0 (non-canonical)
STRAND_CODE = {"1": "+", "2": "-"}  # "0" (undefined) falls through to "."


def opp(strand):
    return {"+": "-", "-": "+", ".": "."}.get(strand, ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sj", required=True, help="STAR {sample}_SJ.out.tab")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--te", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument("--min-unique-reads", type=int, default=1)
    ap.add_argument(
        "--require-canonical", action="store_true",
        help="Only keep junctions with a canonical motif (STAR column 5 != "
        "0). Default here mirrors config chimera.sj_junctions.require_"
        "canonical (true by default, opposite of chimera.reads's).",
    )
    ap.add_argument(
        "--library-strandedness",
        choices=["no", "forward", "reverse", "auto"],
        default="no",
    )
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--te-out", required=False,
        help="Optional second output: only the gene<->TE events (direction "
        "gene_to_te / te_to_gene), same columns as --out.",
    )
    args = ap.parse_args()

    genes = load_bed(args.genes)
    exons = load_bed(args.exons)
    te = load_bed(args.te, n_extra=3)

    tol = max(args.breakpoint_tolerance, 0)
    lib = args.library_strandedness
    if lib == "auto":
        lib = "no"

    gene_meta = {}
    for chrom, (feats, _) in genes.items():
        for s, e, ex in feats:
            gene_meta.setdefault(ex[0], (chrom, s, e, ex[2]))
    te_meta = {}
    for chrom, (feats, _) in te.items():
        for s, e, ex in feats:
            te_meta.setdefault(ex[0], (chrom, s, e, ex[2], ex[3], ex[4], ex[5]))

    def _window(base):
        """Half-open 0-based interval covering the 1-based `base`, +/- tol."""
        return (base - 1 - tol, base + tol)

    rows = []
    n_seen = n_kept = 0

    with open(args.sj) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            n_seen += 1
            chrom = cols[0]
            try:
                intron_start, intron_end = int(cols[1]), int(cols[2])
                motif = int(cols[4])
                unique_reads, multi_reads, overhang = (
                    int(cols[6]), int(cols[7]), int(cols[8])
                )
            except ValueError:
                continue
            annotated = cols[5]
            strand = STRAND_CODE.get(cols[3], ".")

            if unique_reads < args.min_unique_reads:
                continue
            if args.require_canonical and motif not in CANONICAL_MOTIFS:
                continue

            left0, left1 = _window(intron_start - 1)
            right0, right1 = _window(intron_end + 1)

            # Donor = the intron's 5' (transcript-order) splice site,
            # acceptor = its 3'. On "-" strand transcription runs from the
            # higher genomic coordinate to the lower one, so the mapping
            # flips; undefined strand (".") defaults to the "+" mapping
            # (see module docstring -- only reachable when
            # --require-canonical is off).
            if strand == "-":
                donor0, donor1 = right0, right1
                acceptor0, acceptor1 = left0, left1
            else:
                donor0, donor1 = left0, left1
                acceptor0, acceptor1 = right0, right1

            donor_exons = overlapping(exons, chrom, donor0, donor1)
            donor_tes = overlapping(te, chrom, donor0, donor1)
            acceptor_exons = overlapping(exons, chrom, acceptor0, acceptor1)
            acceptor_tes = overlapping(te, chrom, acceptor0, acceptor1)

            donor_genes = sorted({ex[2][0] for ex in donor_exons})
            acceptor_genes = sorted({ex[2][0] for ex in acceptor_exons})
            donor_te_ids = sorted({t[2][0] for t in donor_tes})
            acceptor_te_ids = sorted({t[2][0] for t in acceptor_tes})

            donor_gene_hit = bool(donor_genes)
            donor_te_hit = bool(donor_te_ids)
            acceptor_gene_hit = bool(acceptor_genes)
            acceptor_te_hit = bool(acceptor_te_ids)

            if donor_gene_hit and acceptor_te_hit:
                direction = "gene_to_te"
                gene_id, te_id, gene_side_strand = donor_genes[0], acceptor_te_ids[0], strand
            elif donor_te_hit and acceptor_gene_hit:
                direction = "te_to_gene"
                gene_id, te_id, gene_side_strand = acceptor_genes[0], donor_te_ids[0], strand
            elif donor_gene_hit and acceptor_gene_hit:
                direction = "gene_to_gene"
                gene_id, te_id, gene_side_strand = donor_genes[0], None, strand
            elif donor_te_hit and acceptor_te_hit:
                direction = "te_to_te"
                gene_id, te_id, gene_side_strand = None, donor_te_ids[0], strand
            elif donor_gene_hit:
                direction = "gene_to_other"
                gene_id, te_id, gene_side_strand = donor_genes[0], None, strand
            elif acceptor_gene_hit:
                direction = "other_to_gene"
                gene_id, te_id, gene_side_strand = acceptor_genes[0], None, strand
            elif donor_te_hit:
                direction = "te_to_other"
                gene_id, te_id, gene_side_strand = None, donor_te_ids[0], strand
            elif acceptor_te_hit:
                direction = "other_to_te"
                gene_id, te_id, gene_side_strand = None, acceptor_te_ids[0], strand
            else:
                direction = "other"
                gene_id, te_id, gene_side_strand = None, None, "."

            ambiguous = (donor_gene_hit and donor_te_hit) or (
                acceptor_gene_hit and acceptor_te_hit
            )

            gene_strand = "."
            gene_span = None
            if gene_id is not None and gene_id in gene_meta:
                _, gs, ge, gene_strand = gene_meta[gene_id]
                gene_span = (gs, ge)

            te_family = te_class = te_subfamily = "."
            te_span = None
            if te_id is not None and te_id in te_meta:
                _, ts, tee, _, te_family, te_class, te_subfamily = te_meta[te_id]
                te_span = (ts, tee)

            chimera_type = "."
            antisense = "."
            if direction in ("gene_to_te", "te_to_gene") and gene_span and te_span:
                gs, ge, gst = gene_span[0], gene_span[1], gene_strand
                ts, tee = te_span[0], te_span[1]
                if gst == "+":
                    if tee < gs:
                        chimera_type = "te_initiated"
                    elif ts > ge:
                        chimera_type = "te_terminated"
                    else:
                        chimera_type = "te_exonized"
                elif gst == "-":
                    if ts > ge:
                        chimera_type = "te_initiated"
                    elif tee < gs:
                        chimera_type = "te_terminated"
                    else:
                        chimera_type = "te_exonized"

                te_chrom = chrom  # SJ.out.tab is always cis
                for _s, _e, ex in overlapping(genes, te_chrom, ts, tee):
                    if ex[0] != gene_id and ex[2] == opp(gene_strand):
                        antisense = "yes"
                        break

            transcript_strand = "NA"
            match = "NA"
            if lib == "forward":
                transcript_strand = gene_side_strand
            elif lib == "reverse":
                transcript_strand = opp(gene_side_strand)
            if transcript_strand not in ("NA", ".") and gene_strand != ".":
                match = "yes" if transcript_strand == gene_strand else "no"

            canonical = "yes" if motif in CANONICAL_MOTIFS else "no"
            event_id = f"{chrom}:{intron_start}:{intron_end}:{strand}"
            n_kept += 1

            rows.append([
                event_id, args.sample, chrom, intron_start, intron_end,
                strand, motif, canonical, annotated, unique_reads,
                multi_reads, overhang,
                f"gene:{','.join(donor_genes) or '.'}|te:{','.join(donor_te_ids) or '.'}",
                f"gene:{','.join(acceptor_genes) or '.'}|te:{','.join(acceptor_te_ids) or '.'}",
                direction, "yes" if ambiguous else "no",
                gene_id if gene_id is not None else ".",
                gene_strand,
                te_id if te_id is not None else ".",
                te_subfamily, te_family, te_class, chimera_type, antisense,
                lib, transcript_strand, match,
            ])

    header = [
        "event_id", "sample", "chrom", "intron_start", "intron_end",
        "strand", "motif", "canonical", "annotated", "unique_reads",
        "multi_reads", "overhang", "donor_hits", "acceptor_hits",
        "direction", "direction_ambiguous",
        "gene_id", "gene_strand", "te_id", "te_subfamily", "te_family",
        "te_class", "chimera_type", "antisense_flag", "library_strand",
        "transcript_strand", "gene_strand_match",
    ]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    if args.te_out:
        os.makedirs(os.path.dirname(args.te_out) or ".", exist_ok=True)
        dir_idx = header.index("direction")
        with open_write(args.te_out) as fh:
            fh.write("\t".join(header) + "\n")
            for row in rows:
                if row[dir_idx] in ("gene_to_te", "te_to_gene"):
                    fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{args.sample}: {n_kept}/{n_seen} junctions kept from {args.sj}")


if __name__ == "__main__":
    main()
