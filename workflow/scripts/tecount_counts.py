#!/usr/bin/env python3
"""Merge the per-sample TEcount count tables (TEcount's {sample}.cntTable)
into one feature x sample counts matrix for the sample-QC stage.

The cntTable written by a per-sample TEcount run has a two-column layout:
`gene/TE` (gene ids and TE subfamily names mixed) plus one count column whose
header is the *input BAM path* -- not the sample name -- so samples are mapped
positionally via --sample-names.

Output:

  counts_matrix.tsv  feature x sample integer count matrix (0 where a sample
                     has no reads for the feature). Mirrors the chimera counts
                     matrix layout so the shared sample_qc.R transform/plots
                     modes read it identically.

--feature-class picks which features the matrix keeps (the QC-view scope):
  TE    TE subfamily names only (the default) -- features whose name matches a
        TE GTF element. TEcount builds these keys as `gene_id:family_id:class_id`
        from the TE GTF attributes, so tecount_counts.py reconstructs the same
        keys from --te-gtf.
  gene  gene ids only (from the gene GTF's gene_id attributes).
  all   everything the cntTable holds (genes + TEs, no GTF parsing needed).

The filtering applies only to this QC-view matrix: the per-sample cntTables
are never reduced.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_gtf_keys(path):
    """Feature-key set for one GTF: for the gene GTF that's each feature's
    gene_id; for the TE GTF it's the TEcount element key
    `gene_id:family_id:class_id` (TEcount's TEfeatures.build collapses the
    per-instance transcript ids to these element names)."""
    keys = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = {}
            for name, value in ATTR_RE.findall(fields[8]):
                attrs[name] = value
            if "gene_id" not in attrs:
                continue
            if "family_id" in attrs and "class_id" in attrs:
                # TE GTF line -> TEcount element key
                keys.add(f"{attrs['gene_id']}:{attrs['family_id']}:{attrs['class_id']}")
            else:
                # gene GTF line -> gene_id
                keys.add(attrs["gene_id"])
    return keys


def load_counts(path, feature_keys):
    counts = {}
    with open_read(path) as fh:
        fh.readline()  # header: gene/TE \t <bam path>
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            key = parts[0]
            if feature_keys is not None and key not in feature_keys:
                continue
            counts[key] = int(float(parts[1]))
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--sample-names", required=True, nargs="+")
    ap.add_argument("--gtf", required=False)
    ap.add_argument("--te-gtf", required=False)
    ap.add_argument("--feature-class", required=True,
                    choices=["all", "TE", "gene"])
    ap.add_argument("--out-counts", required=True)
    args = ap.parse_args()

    if len(args.tables) != len(args.sample_names):
        sys.exit("error: --tables and --sample-names must have equal length")

    if args.feature_class == "all":
        feature_keys = None
    elif args.feature_class == "TE":
        if not args.te_gtf:
            sys.exit("error: --feature-class TE requires --te-gtf")
        feature_keys = parse_gtf_keys(args.te_gtf)
    else:
        if not args.gtf:
            sys.exit("error: --feature-class gene requires --gtf")
        feature_keys = parse_gtf_keys(args.gtf)

    features = set()
    sample_counts = []
    for path, sample in zip(args.tables, args.sample_names):
        if not os.path.exists(path):
            sys.exit(f"error: missing count table for sample '{sample}': {path}")
        counts = load_counts(path, feature_keys)
        sample_counts.append((sample, counts))
        features.update(counts)

    with open_write(args.out_counts) as fh:
        fh.write("gene/TE\t" + "\t".join(s for s, _ in sample_counts) + "\n")
        for feature in sorted(features):
            vals = [
                str(counts.get(feature, 0))
                for _, counts in sample_counts
            ]
            fh.write(feature + "\t" + "\t".join(vals) + "\n")

    source = (
        "all"
        if feature_keys is None
        else f"{args.feature_class} ({len(feature_keys)} keys from "
        f"{args.te_gtf if args.feature_class == 'TE' else args.gtf})"
    )
    print(f"{len(features)} {args.feature_class} features across "
          f"{len(args.sample_names)} samples ({source}) -> {args.out_counts}")


if __name__ == "__main__":
    main()
