#!/usr/bin/env python3
"""Classify STAR chimeric junctions as gene-TE chimera candidates, using
read-level breakpoints rather than assembled transcript structure.

The read-evidence counterpart of classify_chimera_assembly.py: same job, same
shape of output, different evidence. Parses STAR's Chimeric.out.junction into
a per-sample, per-event table.

STAR's chimeric detection (--chimOutType Junctions WithinBAM SoftClip, see
rules/align.smk) reports, per chimeric read, one line whose 5' segment is the
"donor" and 3' segment the "acceptor". Each segment's breakpoint is tested
against the gene-exon and TE-insertion tracks (annotation_to_bed.py); a read
whose donor hits a gene exon and acceptor a TE (or vice versa) is a
gene<->TE chimera candidate.

This step is deliberately ANNOTATE-ONLY: every junction is written with its
annotations and flags, and nothing is filtered away here (except, optionally,
events whose STAR junction type is non-canonical when require_canonical_
junction is enabled -- see the config schema; the default keeps everything).
Filtering / evidence decisions are left to the user downstream -- the pipeline
ships the full event table and the counts matrix instead.

Output columns (results/chimera/chimeric_reads/per_sample/{sample}_junctions.tsv.gz):
    event_id, sample, donor_chrom, donor_breakpoint, donor_strand,
    acceptor_chrom, acceptor_breakpoint, acceptor_strand, junction_type,
    canonical, repeat_flag, reads, donor_hits, acceptor_hits, direction,
    direction_ambiguous, gene_id, gene_strand, te_id, te_subfamily,
    te_family, te_class, chimera_type, te_initiated_detail, antisense_flag,
    library_strand, transcript_strand, gene_strand_match, gene_te_distance,
    max_anchor

When --te-out is given, the gene<->TE events (direction gene_to_te /
te_to_gene) are additionally written to that path with the same columns,
so the gene-TE chimeras are available as their own table.

junction_type/canonical: STAR's column-6 value (0 non-canonical .. 6) and a
derived GT/AG-ish yes/no. TE-involved splicing is often non-canonical, so
this is reported but never filtered (see require_canonical config
for the opt-in).

direction: gene_to_te / te_to_gene (the two primary classes -- the gene-TE
chimeras), or one of gene_to_gene / te_to_te / gene_to_other /
other_to_gene / te_to_other / other_to_te / other.  Note that STAR's
"chimeric" is a purely structural call (a read that cannot be explained by
one linear alignment), so the non-gene-TE classes are not just noise: they
also collect circRNA back-splices, read-through transcripts and PCR
chimeras.

direction_ambiguous: "yes" when at least one breakpoint overlaps BOTH a
gene exon and a TE -- TEs sit inside gene bodies routinely, and in that
case several of the direction branches match and the first simply wins.
The reported direction (and the chimera_type derived from it) is then one
defensible reading, not the only one; donor_hits/acceptor_hits carry the
full gene+TE sets for those rows.

chimera_type (gene<->TE events only): BUG FIXED 2026 -- this used to be
decided from where the TE's span sits relative to the GENE's overall
genomic span (TE entirely upstream of the gene -> te_initiated; entirely
downstream -> te_terminated; anywhere else, including squarely inside an
intron -> te_exonized). That is wrong whenever the TE's genomic position
doesn't match the junction's own DIRECTION: a gene_to_te junction into a
TE that happens to sit upstream of the gene's span was scored
te_initiated even though the read shows the gene transcribing INTO the
TE, not the TE initiating anything; a te_to_gene junction from a TE
sitting inside an intron, acting as an alternative promoter that splices
directly into a downstream exon, was scored te_exonized just because the
TE's coordinates fall within the gene's genomic span, even though the
read IS a TE-initiated transcript. This is exactly the bug already fixed
in classify_chimera_splice_junctions.py; see that module's own docstring
for the full reasoning -- both fixes moved to a shared implementation
(chimera_exon_context.py) so the two screens can't drift apart again.

Fixed: chimera_type is now decided from the junction's DIRECTION plus the
gene's own EXON STRUCTURE (from exons.bed, already loaded for donor/
acceptor overlap), matching classify_chimera_assembly.py's vocabulary
exactly (te_initiated / te_terminated / te_exonized -- names unchanged, so
existing consumers of this column keep working), for events where donor
and acceptor are on the SAME chromosome:
  te_to_gene (donor in TE, acceptor in a gene exon) -> te_initiated. The TE
    is transcript-upstream of the exon it splices into, whatever its
    genomic coordinate relative to the gene's overall span.
  gene_to_te (donor in a gene exon, acceptor in TE) -> te_terminated if no
    OTHER annotated exon of the same gene lies transcript-downstream of the
    donor exon (nothing known follows -- the TE plausibly ends the
    transcript); te_exonized if one does (the TE sits inside a region the
    annotation says the gene's transcript continues past -- an internal
    exon, not a true terminus).
For trans events (donor and acceptor on different chromosomes), chimera_type
stays "." -- te_initiated/terminated/exonized would require comparing
coordinates across two different chromosomes.

Rows with direction_ambiguous = "yes" are typed the same way as any other
gene<->TE event above -- direction_ambiguous only flags that the reported
direction was one defensible reading among several (see above), it does not
change how chimera_type is derived from whatever direction was recorded.

te_initiated_detail (te_to_gene events only, "." otherwise): a finer split
that classify_chimera_assembly.py's own te_initiated does NOT distinguish
(so it is reported as its own column, not folded into chimera_type):
  upstream  the acceptor exon is the gene's own most-5' annotated exon (per
            exons.bed) -- the TE splices directly into where the gene
            already starts.
  internal  the gene has an annotated exon further upstream that this
            junction's transcript skips -- the TE is acting as an
            alternative, INTERNAL promoter (e.g. an intronic MT2/MERVL
            case in mouse oocytes/early embryos).

antisense_flag: "yes" when an annotated gene overlaps the TE insertion on the
opposite strand of the assigned gene -- the embedded-TE / sense-antisense
ambiguity class (run twice + IGV curation is the standard treatment for these).

Strand evidence (library_strand / transcript_strand / gene_strand_match):
for stranded libraries the read's aligned strand plus the library type
(forward = read same strand as transcript; reverse = read opposite strand)
gives the transcription strand, which is compared to the annotated gene
strand. For unstranded libraries these columns are NA.

Strand rule (stranded libraries): when a breakpoint overlaps exons of
several genes, a gene on the transcript's own strand is preferred
(prefer_gene_on_strand); a gene<->TE event still on the strand OPPOSITE its
gene is typed antisense_to_gene instead of te_initiated/te_terminated/
te_exonized -- it is antisense transcription through the gene's exon. With
an unstranded library the strand is unknown here (STAR's chimeric strands
are read strands), so no event is re-typed.

gene_te_distance: "trans" (different chromosomes), 0 (TE overlaps the
gene's span) or the gap in bp; "." for events without both a gene and a TE.

max_anchor: for each read, the aligned length of its SHORTER segment (M
bases in the CIGARs, columns 12 and 14); the event reports the best read.
The chimeric-read counterpart of SJ.out.tab's overhang: a breakpoint no read
anchors well on both sides is easier to produce by mis-mapping. On a real
run, local gene-TE events had a median of ~26 bp vs ~18 bp for trans / far
ones. STAR's repeat-length columns (8-9) were checked and are not used: they
measure how far the breakpoint can slide, which tracks the splice motif
rather than artifacts.
"""
import argparse
import bisect
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chimera_exon_context import (
    ANTISENSE_TO_GENE,
    build_gene_exon_positions,
    exon_downstream_of,
    exon_upstream_of,
    prefer_gene_on_strand,
)
from gz_io import open_write


def load_bed(path, n_extra=0):
    """Return {chrom: (feats, max_end)}, where `feats` is a list of
    (start, end, extras_tuple) sorted by start, and `max_end[i]` is the
    running maximum end coordinate over feats[0..i] -- used by overlapping()
    for a correct (not just "usually correct") interval overlap query.
    """
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
        feats.sort()  # by start (tuples compare element-wise)
        running = float("-inf")
        max_end = []
        for _, e, _ in feats:
            running = max(running, e)
            max_end.append(running)
        tracks[chrom] = (feats, max_end)
    return tracks


def overlapping(track, chrom, start0, end0):
    """Features in `track` whose [start, end) overlaps [start0, end0).

    Correct for arbitrarily long/nested/overlapping features, not just the
    common case: `max_end` (built in load_bed) is non-decreasing by
    construction, so bisecting on it finds the leftmost feature whose span
    -- or whose *any preceding* feature's span -- could possibly reach past
    start0, without missing long features that start well before start0 or
    short features sandwiched between them.
    """
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


CANONICAL_TYPES = {1, 2, 3, 4, 5, 6}  # anything but 0 (non-canonical)

_CIGAR_M = re.compile(r"(\d+)M")


def _matched_bases(cigar):
    """Aligned (M) bases in one chimeric segment's CIGAR. STAR writes "-1"
    for an unmapped mate segment, which counts as 0."""
    return sum(int(n) for n in _CIGAR_M.findall(cigar))


def opp(strand):
    return {"+": "-", "-": "+", ".": "."}.get(strand, ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True, help="STAR {sample}_Chimeric.out.junction")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--te", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument(
        "--require-canonical", action="store_true",
        help="Only keep junctions whose STAR junction type is canonical "
        "(GT/AG-ish, type != 0). Default: keep everything and record the type.",
    )
    ap.add_argument(
        "--library-strandedness",
        choices=["no", "forward", "reverse", "auto"],
        default="no",
        help="Per-sample strandedness value as resolved by the workflow "
        "(no = unstranded, forward = read on transcript strand, "
        "reverse = read opposite transcript strand). 'auto' is treated as "
        "unstranded here.",
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

    # Per-gene exon positions, keyed by gene_id -- used below to decide
    # chimera_type from the gene's own exon structure instead of the TE's
    # raw position relative to the gene's overall span. Shared with
    # classify_chimera_splice_junctions.py via chimera_exon_context.py so
    # the two screens' typing logic can't drift apart -- see the module
    # docstring for why the old span-containment test was wrong.
    gene_exon_positions = build_gene_exon_positions(exons)

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

    def donor_locus(bp, strand):
        """Genomic window of the donor segment's LAST ALIGNED base.

        STAR's Chimeric.out.junction column 2 is the first base of the
        DONOR'S INTRON -- not the last base of the aligned segment. The
        aligned base is therefore one step back along the transcript, which
        on a '-' segment means one step forward in genomic coordinates.
        """
        return _window(bp + 1 if strand == "-" else bp - 1)

    def acceptor_locus(bp, strand):
        """Genomic window of the acceptor segment's FIRST ALIGNED base.

        Column 5 is the last base of the ACCEPTOR'S INTRON, so the aligned
        base is one step forward along the transcript (one step back in
        genomic coordinates on a '-' segment).
        """
        return _window(bp - 1 if strand == "-" else bp + 1)

    def read_to_transcript(read_strand):
        """Transcript strand implied by a segment's aligned strand, or "NA"
        when the library is unstranded (STAR's chimeric strands are READ
        strands, unlike SJ.out.tab's motif-derived one)."""
        if lib == "forward":
            return read_strand if read_strand in ("+", "-") else "NA"
        if lib == "reverse":
            return opp(read_strand) if read_strand in ("+", "-") else "NA"
        return "NA"

    def _gene_strand(gene_id):
        return gene_meta.get(gene_id, (None, None, None, "."))[3]

    events = {}

    with open(args.junctions) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            donor_chrom = cols[0]
            acceptor_chrom = cols[3]
            try:
                donor_bp, acceptor_bp = int(cols[1]), int(cols[4])
            except ValueError:
                continue
            donor_strand, acceptor_strand = cols[2], cols[5]
            jtype = cols[6]
            repeat_flag = cols[7]

            d0, d1 = donor_locus(donor_bp, donor_strand)
            a0, a1 = acceptor_locus(acceptor_bp, acceptor_strand)

            donor_exons = overlapping(exons, donor_chrom, d0, d1)
            donor_tes = overlapping(te, donor_chrom, d0, d1)
            acceptor_exons = overlapping(exons, acceptor_chrom, a0, a1)
            acceptor_tes = overlapping(te, acceptor_chrom, a0, a1)

            donor_genes = sorted({ex[2][0] for ex in donor_exons})
            acceptor_genes = sorted({ex[2][0] for ex in acceptor_exons})
            donor_te_ids = sorted({t[2][0] for t in donor_tes})
            acceptor_te_ids = sorted({t[2][0] for t in acceptor_tes})

            donor_gene_hit = bool(donor_genes)
            donor_te_hit = bool(donor_te_ids)
            acceptor_gene_hit = bool(acceptor_genes)
            acceptor_te_hit = bool(acceptor_te_ids)

            # With a stranded library, a gene on the transcript's own strand
            # wins over an overlapping opposite-strand gene (see
            # prefer_gene_on_strand); unstranded keeps the first-sorted gene.
            donor_gene = prefer_gene_on_strand(
                donor_genes, read_to_transcript(donor_strand), _gene_strand)
            acceptor_gene = prefer_gene_on_strand(
                acceptor_genes, read_to_transcript(acceptor_strand), _gene_strand)

            if donor_gene_hit and acceptor_te_hit:
                direction = "gene_to_te"
                gene_id = donor_gene
                te_id = acceptor_te_ids[0]
                gene_side_strand = donor_strand
            elif donor_te_hit and acceptor_gene_hit:
                direction = "te_to_gene"
                gene_id = acceptor_gene
                te_id = donor_te_ids[0]
                gene_side_strand = acceptor_strand
            elif donor_gene_hit and acceptor_gene_hit:
                direction = "gene_to_gene"
                gene_id = donor_gene
                te_id = None
                gene_side_strand = donor_strand
            elif donor_te_hit and acceptor_te_hit:
                direction = "te_to_te"
                gene_id = None
                te_id = donor_te_ids[0]
                gene_side_strand = donor_strand
            elif donor_gene_hit:
                direction = "gene_to_other"
                gene_id = donor_gene
                te_id = None
                gene_side_strand = donor_strand
            elif acceptor_gene_hit:
                direction = "other_to_gene"
                gene_id = acceptor_gene
                te_id = None
                gene_side_strand = acceptor_strand
            elif donor_te_hit:
                direction = "te_to_other"
                gene_id = None
                te_id = donor_te_ids[0]
                gene_side_strand = donor_strand
            elif acceptor_te_hit:
                direction = "other_to_te"
                gene_id = None
                te_id = acceptor_te_ids[0]
                gene_side_strand = acceptor_strand
            else:
                direction = "other"
                gene_id = None
                te_id = None
                gene_side_strand = "."

            # A single breakpoint can overlap BOTH a gene exon and a TE
            # (TEs sit inside gene bodies all the time), in which case more
            # than one branch above would have matched and the first one
            # simply won.  Flag those rows: `direction` (and the
            # chimera_type derived from it) is then one defensible reading
            # of an ambiguous locus, not the only one -- check
            # donor_hits/acceptor_hits, which carry the full gene+TE sets.
            ambiguous = (donor_gene_hit and donor_te_hit) or (
                acceptor_gene_hit and acceptor_te_hit
            )

            key = (donor_chrom, donor_bp, donor_strand, acceptor_chrom,
                   acceptor_bp, acceptor_strand, direction)
            ev = events.setdefault(key, {"reads": 0, "gene_id": gene_id,
                                          "te_id": te_id, "max_anchor": 0})
            ev["reads"] += 1
            # the read's shorter segment, in aligned bases (CIGAR columns
            # 12 / 14); max over the event's reads -- see max_anchor in the
            # module docstring
            if len(cols) >= 14:
                ev["max_anchor"] = max(ev["max_anchor"], min(
                    _matched_bases(cols[11]), _matched_bases(cols[13])))
            if ev["reads"] == 1:
                ev.update(
                    {
                        "donor_hits": ",".join(donor_genes) or ".",
                        "acceptor_hits": ",".join(acceptor_genes) or ".",
                        "donor_te_hits": ",".join(donor_te_ids) or ".",
                        "acceptor_te_hits": ",".join(acceptor_te_ids) or ".",
                        "junction_type": jtype,
                        "repeat_flag": repeat_flag,
                        "gene_side_strand": gene_side_strand,
                        "ambiguous": "yes" if ambiguous else "no",
                    }
                )
            else:
                # first-seen gene/TE win for stability; keep the union of hits
                ev["donor_hits"] = _union(ev["donor_hits"], donor_genes)
                ev["acceptor_hits"] = _union(ev["acceptor_hits"], acceptor_genes)
                ev["donor_te_hits"] = _union(ev["donor_te_hits"], donor_te_ids)
                ev["acceptor_te_hits"] = _union(ev["acceptor_te_hits"], acceptor_te_ids)

    # Resolve per-event annotations (gene span/strand, TE family/class, type,
    # antisense, strand evidence).
    rows = []
    for (donor_chrom, donor_bp, donor_strand, acceptor_chrom, acceptor_bp,
         acceptor_strand, direction), ev in sorted(events.items()):
        if args.require_canonical:
            try:
                if int(ev["junction_type"]) not in CANONICAL_TYPES:
                    continue
            except ValueError:
                continue
        gene_id = ev["gene_id"]
        te_id = ev["te_id"]

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
        te_initiated_detail = "."
        antisense = "."
        same_chrom = donor_chrom == acceptor_chrom
        # Trans events (different chromosomes) leave chimera_type as "." --
        # te_initiated/terminated/exonized would require comparing
        # coordinates across two different chromosomes.
        if direction in ("gene_to_te", "te_to_gene") and gene_span and te_span and same_chrom:
            gst = gene_strand

            if direction == "te_to_gene":
                # Donor in TE, acceptor in a gene exon: the TE is
                # transcript-upstream of the exon it splices into, whatever
                # its genomic coordinate relative to the gene's overall
                # span -- see the module docstring for why span-containment
                # was wrong here.
                chimera_type = "te_initiated"
                a0, a1 = acceptor_locus(acceptor_bp, acceptor_strand)
                accept_pos = (a0 + a1) // 2
                te_initiated_detail = (
                    "internal"
                    if exon_upstream_of(gene_exon_positions, gene_id, accept_pos, gst)
                    else "upstream"
                )
            else:  # gene_to_te
                # Donor in a gene exon, acceptor in TE: terminated only if
                # no OTHER annotated exon of this gene lies further
                # downstream than the donor exon (nothing known follows the
                # TE); exonized if one does (the annotation says the gene's
                # transcript continues past this point).
                d0, d1 = donor_locus(donor_bp, donor_strand)
                donor_pos = (d0 + d1) // 2
                chimera_type = (
                    "te_exonized"
                    if exon_downstream_of(gene_exon_positions, gene_id, donor_pos, gst)
                    else "te_terminated"
                )

        if direction in ("gene_to_te", "te_to_gene") and gene_span and te_span:
            # antisense: annotated gene overlapping the TE insertion on the
            # strand opposite the assigned gene.  TE is on the acceptor side
            # for gene_to_te, donor side for te_to_gene -- te_chrom is
            # always correct regardless of same_chrom, so this still runs
            # for trans events too.
            te_chrom = acceptor_chrom if direction == "gene_to_te" else donor_chrom
            for _s, _e, ex in overlapping(genes, te_chrom, te_span[0], te_span[1]):
                if ex[0] != gene_id and ex[2] == opp(gene_strand):
                    antisense = "yes"
                    break

        # strand evidence
        transcript_strand = read_to_transcript(ev["gene_side_strand"])
        match = "NA"
        if transcript_strand != "NA" and gene_strand in ("+", "-"):
            match = "yes" if transcript_strand == gene_strand else "no"

        # Strand rule (stranded libraries only -- unstranded leaves match
        # NA and the type unchanged): a transcript on the strand opposite
        # the assigned gene is antisense transcription through that gene's
        # exon, not initiation/termination/exonization of the gene.
        if match == "no" and chimera_type in (
                "te_initiated", "te_terminated", "te_exonized"):
            chimera_type = ANTISENSE_TO_GENE
            te_initiated_detail = "."

        # Genomic distance between the gene and the TE: "trans" on
        # different chromosomes, 0 when the TE overlaps the gene's span,
        # otherwise the gap in bp. Most chimeric-read gene<->TE events on a
        # real run joined a gene to a TE on another chromosome or >200 kb
        # away -- the partner pattern of template switching / chimeric
        # ligation -- so this lets them be separated from local events.
        gene_te_distance = "."
        if gene_span and te_span:
            gene_chrom = gene_meta[gene_id][0]
            te_chrom_ = te_meta[te_id][0]
            if gene_chrom != te_chrom_:
                gene_te_distance = "trans"
            else:
                gene_te_distance = max(0, te_span[0] - gene_span[1],
                                       gene_span[0] - te_span[1])

        try:
            canonical = "yes" if int(ev["junction_type"]) in CANONICAL_TYPES else "no"
        except ValueError:
            canonical = "no"
        event_id = f"{donor_chrom}:{donor_bp}:{donor_strand}:{acceptor_chrom}:{acceptor_bp}:{acceptor_strand}:{direction}"

        rows.append(
            [
                event_id, args.sample, donor_chrom, donor_bp, donor_strand,
                acceptor_chrom, acceptor_bp, acceptor_strand,
                ev["junction_type"], canonical,
                ev["repeat_flag"], ev["reads"],
                f"gene:{ev['donor_hits']}|te:{ev['donor_te_hits']}",
                f"gene:{ev['acceptor_hits']}|te:{ev['acceptor_te_hits']}",
                direction,
                ev["ambiguous"],
                gene_id if gene_id is not None else ".",
                gene_strand,
                te_id if te_id is not None else ".",
                te_subfamily, te_family, te_class, chimera_type,
                te_initiated_detail, antisense,
                lib, transcript_strand, match, gene_te_distance,
                ev["max_anchor"],
            ]
        )

    header = [
        "event_id", "sample", "donor_chrom", "donor_breakpoint", "donor_strand",
        "acceptor_chrom", "acceptor_breakpoint", "acceptor_strand",
        "junction_type", "canonical",
        "repeat_flag", "reads", "donor_hits", "acceptor_hits", "direction",
        "direction_ambiguous",
        "gene_id", "gene_strand", "te_id", "te_subfamily", "te_family",
        "te_class",
        "chimera_type", "te_initiated_detail", "antisense_flag",
        "library_strand", "transcript_strand",
        "gene_strand_match", "gene_te_distance", "max_anchor",
    ]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    if args.te_out:
        os.makedirs(os.path.dirname(args.te_out), exist_ok=True)
        dir_idx = header.index("direction")
        with open_write(args.te_out) as fh:
            fh.write("\t".join(header) + "\n")
            for row in rows:
                if row[dir_idx] in ("gene_to_te", "te_to_gene"):
                    fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{args.sample}: {len(rows)} events from {args.junctions}")


def _union(csv, ids):
    have = set(x for x in csv.split(",") if x and x != ".")
    have.update(ids)
    return ",".join(sorted(have))


if __name__ == "__main__":
    main()
