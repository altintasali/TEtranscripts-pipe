"""Shared exon-neighbor logic for deciding chimera_type from junction
DIRECTION plus a gene's own annotated exon structure, used identically by
classify_chimera_chimeric_reads.py and classify_chimera_splice_junctions.py
so the two screens cannot drift apart on what te_initiated / te_terminated /
te_exonized mean. Both used to decide chimera_type from where the TE's span
sits relative to the GENE's overall genomic span, which is wrong for a TE
acting as an intronic alternative promoter or terminator -- see either
classifier's own module docstring for the full reasoning.
"""


def build_gene_exon_positions(exons_track):
    """{gene_id: [(start, end), ...]} from a load_bed()-style exons track
    (exons.bed is gene_id-keyed -- see annotation_to_bed.py: a flattened
    union of every transcript's exons under that gene)."""
    gene_exon_positions = {}
    for _chrom, (feats, _max_end) in exons_track.items():
        for s, e, ex in feats:
            gene_exon_positions.setdefault(ex[0], []).append((s, e))
    return gene_exon_positions


def exon_upstream_of(gene_exon_positions, gene_id, pos, gene_strand):
    """True if gene_id has an annotated exon entirely transcript-upstream
    (5') of genomic position pos, strand-aware."""
    for s, e in gene_exon_positions.get(gene_id, []):
        if gene_strand == "-":
            if s > pos:
                return True
        elif e < pos:
            return True
    return False


def exon_downstream_of(gene_exon_positions, gene_id, pos, gene_strand):
    """True if gene_id has an annotated exon entirely transcript-
    downstream (3') of genomic position pos, strand-aware."""
    for s, e in gene_exon_positions.get(gene_id, []):
        if gene_strand == "-":
            if e < pos:
                return True
        elif s > pos:
            return True
    return False


# chimera_type for a gene<->TE call whose transcript runs on the strand
# OPPOSITE the assigned gene. Such a call is antisense transcription through
# the gene's exon, not initiation/termination/exonization of that gene, so
# it gets its own label instead of te_initiated/te_terminated/te_exonized
# (which describe the GENE's transcript). Shared by all three classifiers.
#
# Not to be confused with the older antisense_flag column: that one says
# ANOTHER annotated gene overlaps the TE on the opposite strand, and says
# nothing about which strand the transcript itself was on.
ANTISENSE_TO_GENE = "antisense_to_gene"


def prefer_gene_on_strand(gene_ids, transcript_strand, strand_of):
    """Pick the gene for a breakpoint that overlaps exons of several genes.

    gene_ids is the sorted candidate list the classifiers already build;
    strand_of(gene_id) returns that gene's strand. When the transcript's
    strand is known ("+"/"-"), the first gene ON that strand wins -- a
    breakpoint inside two overlapping genes on opposite strands used to go
    to whichever sorted first, which could turn an ordinary sense splice
    into an apparent antisense one. With no same-strand candidate (or an
    unknown transcript strand) the old first-sorted choice is kept, and the
    caller decides whether the result is antisense_to_gene.
    """
    if not gene_ids:
        return None
    if transcript_strand in ("+", "-"):
        for gid in gene_ids:
            if strand_of(gid) == transcript_strand:
                return gid
    return gene_ids[0]


# Known gene structure that merely touches a TE -- kept as its own class so
# it stays visible, but not a gene-TE chimera call:
#   annotated_splice                      an SJ junction that is an annotated
#                                         intron of the reference GTF
#                                         (classify_chimera_splice_junctions.py)
#   annotated_terminal_exon_embedded_te   an assembled transcript whose
#                                         TE-overlapping last exon is an
#                                         annotated last exon, i.e. a TE in an
#                                         ordinary 3' UTR -- the 3' counterpart
#                                         of annotated_promoter_embedded_te
#                                         (classify_chimera_assembly.py)
ANNOTATED_SPLICE = "annotated_splice"
ANNOTATED_TERMINAL_EXON = "annotated_terminal_exon_embedded_te"
