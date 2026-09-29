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
