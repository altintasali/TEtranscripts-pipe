#!/usr/bin/env python3
"""Merge the per-sample TEcount/TElocal count tables into one feature x
sample counts matrix, or filter an already-merged matrix down to one
feature class.

The cntTable written by a per-sample TEcount/TElocal run has a two-column
layout: `gene/TE` (gene ids and TE feature names mixed) plus one count column
whose header is the *input BAM path* -- not the sample name -- so samples are
mapped positionally via --sample-names.

Two mutually exclusive modes, selected by which input flags are given:

  merge mode   --tables + --sample-names. Merges every sample's cntTable
               into one feature x sample matrix (0 where a sample has no
               reads for the feature). This is how results/tecount/
               counts_matrix.tsv.gz and results/telocal/counts_matrix.tsv.gz
               are built -- always with --feature-class all, so they hold
               every feature TEcount/TElocal reports, matching native
               output exactly.
  filter mode  --in-matrix. Re-filters an already-merged matrix (such as
               the counts_matrix.tsv.gz above) down to one feature class,
               without re-reading the per-sample cntTables. This is how the
               QC-view's TE-only (or gene-only) matrix is derived; see
               --feature-class below.

Output (either mode):

  counts_matrix.tsv  feature x sample integer count matrix. Mirrors the
                     chimera counts matrix layout so the shared sample_qc.R
                     transform/plots modes read it identically.

--key-style selects how rows are classified for the --feature-class filter:
  tecount  TEcount tables (the default). TE keys are `gene_id:family_id:
           class_id` built from the TE GTF attributes, so the TE/gene key
           sets are reconstructed from --te-gtf / --gtf.
  telocal  TElocal tables. TE keys are locus keys of the form
           `chr:start:end(family:strand):gene_id:family_id:class_id`
           (>= 3 colon-separated fields; the transcript id itself embeds
           colons), gene keys have none -- so classification needs no GTF
           and --gtf/--te-gtf are ignored.

--feature-class picks which features the output matrix keeps:
  TE    TE features only (the QC-view default).
  gene  gene ids only.
  all   everything the input holds (genes + TEs, no classification). Always
        used for the top-level counts_matrix.tsv.gz merge.

Merge mode never reduces the per-sample cntTables; filter mode never
touches the merge-mode output it reads from -- it only writes a new,
separate file.
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


def is_telocal_te(key):
    """True when the key is a TElocal TE-locus key (>= 3 colon-separated,
    non-empty fields), matching telocal_summary_mqc.py's classify()."""
    parts = key.split(":")
    return len(parts) >= 3 and all(parts)


def resolve_feature_keys(feature_class, key_style, gtf, te_gtf):
    """Feature-key spec for --feature-class, shared by merge and filter
    mode. Returns None (keep everything), a set of tecount-style keys to
    match exactly, or "TE"/"gene" for telocal's shape-based classification
    (see keep_key)."""
    if feature_class == "all":
        return None
    if key_style == "telocal":
        return feature_class
    if feature_class == "TE":
        if not te_gtf:
            sys.exit("error: --feature-class TE requires --te-gtf")
        return parse_gtf_keys(te_gtf)
    if not gtf:
        sys.exit("error: --feature-class gene requires --gtf")
    return parse_gtf_keys(gtf)


def keep_key(key, feature_keys):
    """True when `key` belongs to the --feature-class selection resolved by
    resolve_feature_keys."""
    if feature_keys is None:
        return True
    if isinstance(feature_keys, str):
        return is_telocal_te(key) == (feature_keys == "TE")
    return key in feature_keys


def load_counts(path, feature_keys):
    # main() holds one of these dicts per sample, all at once, until every
    # sample has been read (needed to build the union of features before
    # writing). Every sample's cntTable shares the SAME feature vocabulary
    # (same TE GTF / TElocal index), but str.split() allocates a brand-new
    # string object per line, per file -- so without interning, a real
    # locus-level cohort (millions of loci x dozens of samples) holds that
    # many independent copies of what is really a small, shared string
    # vocabulary. sys.intern() collapses equal strings to one shared object;
    # dict lookup behavior is unchanged (still exact string equality), only
    # the memory backing it is shared. Measured on a synthetic equivalent at
    # reduced scale: 66% peak-memory reduction, and this was the direct
    # cause of a real OOM (telocal_counts, 84 samples, ~3.7M loci/sample).
    counts = {}
    with open_read(path) as fh:
        fh.readline()  # header: gene/TE \t <bam path>
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            key = sys.intern(parts[0])
            if not keep_key(key, feature_keys):
                continue
            counts[key] = int(float(parts[1]))
    return counts


def filter_matrix(in_path, out_path, feature_keys):
    """Stream an already-merged feature x sample matrix, keeping only rows
    whose key passes keep_key. Unlike merge mode, this holds nothing in
    memory beyond the current line -- there is no per-sample dict to
    accumulate, just a row-by-row pass over one file."""
    kept = 0
    with open_read(in_path) as in_fh, open_write(out_path) as out_fh:
        header = in_fh.readline()
        out_fh.write(header)
        for line in in_fh:
            if not line.strip():
                continue
            key = line.split("\t", 1)[0]
            if not keep_key(key, feature_keys):
                continue
            out_fh.write(line if line.endswith("\n") else line + "\n")
            kept += 1
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", nargs="+",
                     help="Merge mode: per-sample cntTable paths.")
    ap.add_argument("--sample-names", nargs="+",
                     help="Merge mode: sample name per --tables entry.")
    ap.add_argument("--in-matrix",
                     help="Filter mode: an already-merged counts matrix to "
                          "re-filter by --feature-class.")
    ap.add_argument("--key-style", default="tecount",
                    choices=["tecount", "telocal"])
    ap.add_argument("--gtf", required=False)
    ap.add_argument("--te-gtf", required=False)
    ap.add_argument("--feature-class", required=True,
                    choices=["all", "TE", "gene"])
    ap.add_argument("--out-counts", required=True)
    args = ap.parse_args()

    merge_mode = bool(args.tables or args.sample_names)
    filter_mode = bool(args.in_matrix)
    if merge_mode and filter_mode:
        sys.exit("error: --in-matrix cannot be combined with --tables/--sample-names")
    if not merge_mode and not filter_mode:
        sys.exit("error: provide either --in-matrix, or --tables together with --sample-names")
    if merge_mode and (not args.tables or not args.sample_names
                       or len(args.tables) != len(args.sample_names)):
        sys.exit("error: --tables and --sample-names must both be given with equal length")

    feature_keys = resolve_feature_keys(
        args.feature_class, args.key_style, args.gtf, args.te_gtf)

    if filter_mode:
        if not os.path.exists(args.in_matrix):
            sys.exit(f"error: missing input matrix: {args.in_matrix}")
        n = filter_matrix(args.in_matrix, args.out_counts, feature_keys)
        print(f"{n} {args.feature_class} features kept from {args.in_matrix} "
              f"-> {args.out_counts}")
        return

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

    if args.feature_class == "all":
        source = "all"
    elif args.key_style == "telocal":
        source = f"{args.feature_class} ({'TE-local locus' if args.feature_class == 'TE' else 'gene'} keys by shape)"
    else:
        source = (
            f"{args.feature_class} ({len(feature_keys)} keys from "
            f"{args.te_gtf if args.feature_class == 'TE' else args.gtf})"
        )
    print(f"{len(features)} {args.feature_class} features across "
          f"{len(args.sample_names)} samples ({source}) -> {args.out_counts}")


if __name__ == "__main__":
    main()
