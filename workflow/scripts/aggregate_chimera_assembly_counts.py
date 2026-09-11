#!/usr/bin/env python3
"""Sum results/chimera/assembly/counts_matrix.tsv.gz (transcript_id x sample,
StringTie-derived raw count ESTIMATES, see quantify_chimera_assembly.py) down
to one row per (matched_gene_id, te_id, chimera_type) -- the grain DESeq2/
edgeR should actually test at, not the raw MSTRG transcript_id.

Why this grain and not transcript_id or (gene, te) alone:

  * transcript_id is too fragmented for a DE test. StringTie routinely
    assembles one real transcriptional event as several near-identical
    isoforms differing by a few bp at an imprecise exon boundary (see
    classify_chimera_assembly.py's breakpoint-tolerance discussion) --
    splitting the read support for one event across N noisy rows.
  * Summing every isoform of a (gene_id, te_id) pair together, ignoring
    chimera_type, is wrong: a single candidate pair can carry isoforms of
    DIFFERENT chimera_type (te_initiated vs te_terminated vs te_exonized),
    which are mechanistically distinct regulatory events -- pooling them
    can cancel out or mask a real signal specific to just one mechanism.
  * (gene_id, te_id, chimera_type) only pools the sample dimension (never
    collapsed) and only pools isoforms that already share the same
    structural class -- the legitimate unit.

Rows whose transcripts.tsv.gz candidate has matched_gene_id=="." and/or
te_id=="." (te_initiated_intergenic has no matched gene; unspliced_te_only
has no splice evidence) are excluded, same filter idiom chimera_evidence.py
uses for its own (gene, te) pair aggregation.

Two outputs share the gene_te_chimera_id row key:

  --out-counts (gene_te_chimera_id x sample, DESeq2/edgeR-ready):
    gene_te_chimera_id   f"{matched_gene_id}:{te_id}:{chimera_type}" -- a
                         single joined key column (not three separate
                         columns) so read.delim(file, row.names=1) in R
                         hands DESeqDataSetFromMatrix a pure numeric matrix
                         with no further munging.
    {sample columns}     summed raw counts (integers), one column per
                         sample, in counts_matrix.tsv.gz's own header order.

  --out-annotation (one row per gene_te_chimera_id -- the rowData a DESeq2
    user attaches to label result tables/volcano plots, since counts_matrix
    alone carries no gene symbol, TE class/family/subfamily, or locus):
    gene_te_chimera_id, gene_id, gene_symbol, gene_locus, te_id, te_locus,
    te_subfamily, te_family, te_class, chimera_type, n_transcripts,
    assembly_transcript_ids.
    te_subfamily/te_family/te_class come from the curated TE GTF (joined in
    classify_chimera_assembly.py), not from StringTie, so they are
    IDENTICAL across every transcript sharing a te_id -- taking them from
    any one member transcript is exact, not an approximation.
    gene_symbol falls back to gene_id when --gene-names is omitted or the
    id is unresolved, same convention as chimera_candidates_table_mqc.py /
    chimera_candidates_explorer.R. gene_locus/te_locus are ready-to-paste
    IGV coordinates ("chr:start-end", 1-based inclusive) joined from
    genes.bed/te.bed -- same BED-to-locus convention (and off-by-one fix:
    BED start is 0-based half-open) as
    chimera_candidates_explorer.R's read_bed_loci(); "." when the id isn't
    in the BED.

Usage:
  aggregate_chimera_assembly_counts.py \\
      --transcripts results/chimera/assembly/transcripts.tsv.gz \\
      --counts results/chimera/assembly/counts_matrix.tsv.gz \\
      --gene-names results/reference/gene_id_to_name.tsv.gz \\
      --genes-bed results/reference/genes.bed \\
      --te-bed results/reference/te.bed \\
      --out-counts gene_te_chimera_counts_matrix.tsv.gz \\
      --out-annotation gene_te_chimera_counts_annotation.tsv.gz
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def load_gene_symbols(path):
    """{gene_id: gene_symbol}, empty when --gene-names is omitted or the
    file doesn't exist -- same fallback convention as
    chimera_candidates_table_mqc.py's load_symbols()."""
    symbols = {}
    if not path or not os.path.exists(path):
        return symbols
    with open_read(path) as fh:
        fh.readline()
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1] not in (".", ""):
                symbols[fields[0]] = fields[1]
    return symbols


def load_bed_loci(path):
    """{id: "chrom:start-end"} from a BED file keyed on column 4 (chrom,
    start, end, id, score, strand, ...; no header row -- see
    annotation_to_bed.py). BED start is 0-based half-open; IGV's locus box
    is 1-based inclusive, hence the +1 on start only -- same convention as
    chimera_candidates_explorer.R's read_bed_loci()."""
    loci = {}
    with open_read(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, start, end, id_ = cols[0], int(cols[1]), cols[2], cols[3]
            loci[id_] = f"{chrom}:{start + 1}-{end}"
    return loci


def load_transcript_groups(path):
    """{transcript_id: None or (gene, te, chimera_type, te_subfamily,
    te_family, te_class)} -- None for rows that are not real gene-TE
    candidates (matched_gene_id or te_id == "."), same exclusion
    chimera_evidence.py applies. Every transcript_id is still a dict key
    (even when its value is None) so the counts-matrix pass below can tell
    "known, filtered out" apart from "unknown id", which would mean the two
    input files are out of sync."""
    groups = {}
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            tid = cols[idx["transcript_id"]]
            gene = cols[idx["matched_gene_id"]]
            te = cols[idx["te_id"]]
            ctype = cols[idx["chimera_type"]]
            if gene in (".", "") or te in (".", "") or ctype in (".", ""):
                groups[tid] = None
            else:
                groups[tid] = (
                    gene, te, ctype,
                    cols[idx["te_subfamily"]], cols[idx["te_family"]],
                    cols[idx["te_class"]],
                )
    return groups


def aggregate(counts_path, groups):
    """{(gene, te, chimera_type): [per-sample summed int]} and
    {(gene, te, chimera_type): [transcript_id, ...]} (the members summed
    into that group -- doubles as n_transcripts and the annotation table's
    assembly_transcript_ids), plus the sample list (from counts_matrix.tsv.gz's
    own header) and include/exclude counts for the summary line."""
    sums = {}
    members = {}
    n_included = n_excluded = 0
    with open_read(counts_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            tid = parts[0]
            if tid not in groups:
                sys.exit(f"error: transcript_id {tid!r} in {counts_path} not "
                         f"found in --transcripts -- inputs out of sync")
            info = groups[tid]
            if info is None:
                n_excluded += 1
                continue
            key = info[:3]
            n_included += 1
            values = [int(v) for v in parts[1:]]
            row = sums.setdefault(key, [0] * len(samples))
            for i, v in enumerate(values):
                row[i] += v
            members.setdefault(key, []).append(tid)
    return samples, sums, members, n_included, n_excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcripts", required=True,
                     help="results/chimera/assembly/transcripts.tsv.gz")
    ap.add_argument("--counts", required=True,
                     help="results/chimera/assembly/counts_matrix.tsv.gz")
    ap.add_argument("--gene-names", default=None,
                     help="results/reference/gene_id_to_name.tsv.gz, to "
                     "label the annotation table's gene_symbol column")
    ap.add_argument("--genes-bed", required=True,
                     help="results/reference/genes.bed")
    ap.add_argument("--te-bed", required=True,
                     help="results/reference/te.bed")
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-annotation", required=True)
    args = ap.parse_args()

    groups = load_transcript_groups(args.transcripts)
    symbols = load_gene_symbols(args.gene_names)
    gene_loci = load_bed_loci(args.genes_bed)
    te_loci = load_bed_loci(args.te_bed)
    samples, sums, members, n_included, n_excluded = aggregate(args.counts, groups)
    # te_subfamily/family/class for a group: read them off any one member
    # transcript (guaranteed identical across the group -- see
    # load_transcript_groups's docstring) rather than threading a third
    # dict through aggregate().
    te_annot = {key: groups[tids[0]][3:] for key, tids in members.items()}

    os.makedirs(os.path.dirname(os.path.abspath(args.out_counts)), exist_ok=True)
    with open_write(args.out_counts) as fh:
        fh.write("gene_te_chimera_id\t" + "\t".join(samples) + "\n")
        for gene, te, ctype in sorted(sums):
            row = sums[(gene, te, ctype)]
            fh.write(f"{gene}:{te}:{ctype}\t" + "\t".join(str(v) for v in row) + "\n")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_annotation)), exist_ok=True)
    with open_write(args.out_annotation) as fh:
        fh.write("gene_te_chimera_id\tgene_id\tgene_symbol\tgene_locus\tte_id\t"
                 "te_locus\tte_subfamily\tte_family\tte_class\tchimera_type\t"
                 "n_transcripts\tassembly_transcript_ids\n")
        for gene, te, ctype in sorted(sums):
            key = (gene, te, ctype)
            tids = sorted(members[key])
            te_subfamily, te_family, te_class = te_annot[key]
            gene_symbol = symbols.get(gene, gene)
            fh.write("\t".join([
                f"{gene}:{te}:{ctype}", gene, gene_symbol,
                gene_loci.get(gene, "."), te, te_loci.get(te, "."),
                te_subfamily, te_family, te_class, ctype,
                str(len(tids)), ",".join(tids),
            ]) + "\n")

    print(f"{len(sums)} (gene,te,chimera_type) groups from {n_included} "
          f"candidate transcripts ({n_excluded} excluded: no matched gene "
          f"and/or TE) -> {args.out_counts}, {args.out_annotation}")


if __name__ == "__main__":
    main()
