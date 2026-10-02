#!/usr/bin/env python3
"""Sum results/chimera/splice_junctions/counts_matrix.tsv.gz (event_id x
sample, STAR's unique_reads per junction, see
chimera_splice_junctions_counts.py) down to one row per (gene_id, te_id,
chimera_type) -- the same grain and the same row key as the assembly
screen's gene_te_chimera_counts_matrix.tsv.gz
(aggregate_chimera_assembly_counts.py, whose docstring explains why this
grain), so a hit from one matrix can be looked up in the other.

The assembly matrix is the one to test with; this one is for confirming a
hit with exact junction reads. Do not sum the two: they count overlapping
reads.

Caveat: a read crossing two gene-TE junctions of the same pair (gene exon
-> TE exon -> gene exon, i.e. te_exonized) is counted once per junction.

Junction events whose gene_id, te_id or chimera_type is "." are left out,
same as the assembly aggregation.

Two outputs share the gene_te_chimera_id row key:

  --out-counts (gene_te_chimera_id x sample, DESeq2/edgeR-ready):
    gene_te_chimera_id   f"{gene_id}:{te_id}:{chimera_type}"
    {sample columns}     summed unique junction reads, one column per
                         sample, in counts_matrix.tsv.gz's own header order.

  --out-annotation (one row per gene_te_chimera_id):
    gene_te_chimera_id, gene_id, gene_symbol, gene_locus, te_id, te_locus,
    te_subfamily, te_family, te_class, chimera_type, n_junctions,
    sj_event_ids.

Usage:
  aggregate_chimera_splice_junctions_counts.py \\
      --junctions results/chimera/splice_junctions/te-gene-junctions.tsv.gz \\
      --counts results/chimera/splice_junctions/counts_matrix.tsv.gz \\
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
from aggregate_chimera_assembly_counts import load_bed_loci, load_gene_symbols
from gz_io import open_read, open_write


def load_junction_groups(path):
    """{event_id: None or (gene, te, chimera_type, te_subfamily, te_family,
    te_class)} -- None for events that are not typed gene-TE junctions.
    Every event_id is a key, so an id in the counts matrix but not here
    means the two inputs are out of sync."""
    groups = {}
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            eid = cols[idx["event_id"]]
            gene = cols[idx["gene_id"]]
            te = cols[idx["te_id"]]
            ctype = cols[idx["chimera_type"]]
            if gene in (".", "") or te in (".", "") or ctype in (".", ""):
                groups[eid] = None
            else:
                groups[eid] = (
                    gene, te, ctype,
                    cols[idx["te_subfamily"]], cols[idx["te_family"]],
                    cols[idx["te_class"]],
                )
    return groups


def aggregate(counts_path, groups):
    """Per-sample sums and member event_ids per (gene, te, chimera_type),
    plus the sample list and include/exclude counts."""
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
            eid = parts[0]
            if eid not in groups:
                sys.exit(f"error: event_id {eid!r} in {counts_path} not "
                         f"found in --junctions -- inputs out of sync")
            info = groups[eid]
            if info is None:
                n_excluded += 1
                continue
            key = info[:3]
            n_included += 1
            row = sums.setdefault(key, [0] * len(samples))
            for i, v in enumerate(parts[1:]):
                row[i] += int(v)
            members.setdefault(key, []).append(eid)
    return samples, sums, members, n_included, n_excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True,
                     help="results/chimera/splice_junctions/te-gene-junctions.tsv.gz")
    ap.add_argument("--counts", required=True,
                     help="results/chimera/splice_junctions/counts_matrix.tsv.gz")
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

    groups = load_junction_groups(args.junctions)
    symbols = load_gene_symbols(args.gene_names)
    gene_loci = load_bed_loci(args.genes_bed)
    te_loci = load_bed_loci(args.te_bed)
    samples, sums, members, n_included, n_excluded = aggregate(args.counts, groups)
    te_annot = {key: groups[eids[0]][3:] for key, eids in members.items()}

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
                 "n_junctions\tsj_event_ids\n")
        for gene, te, ctype in sorted(sums):
            key = (gene, te, ctype)
            eids = sorted(members[key])
            te_subfamily, te_family, te_class = te_annot[key]
            fh.write("\t".join([
                f"{gene}:{te}:{ctype}", gene, symbols.get(gene, gene),
                gene_loci.get(gene, "."), te, te_loci.get(te, "."),
                te_subfamily, te_family, te_class, ctype,
                str(len(eids)), ",".join(eids),
            ]) + "\n")

    print(f"{len(sums)} (gene,te,chimera_type) groups from {n_included} "
          f"junction events ({n_excluded} excluded: no gene, TE or type) "
          f"-> {args.out_counts}, {args.out_annotation}")


if __name__ == "__main__":
    main()
