#!/usr/bin/env python3
"""Build transcript_id x sample TPM/counts/CPM matrices for chimera-assembly
candidates, pulling values out of each sample's -e/-B re-quantified StringTie
GTF (results/chimera/assembly/per_sample/quant/{sample}.transcripts.gtf, see
chimera_assembly.smk's stringtie_requantify rule).

Kept separate from classify_chimera_assembly.py the same way chimera_reads_counts.py
is separate from classify_chimera_reads.py: structural classification and
expression aggregation are different concerns, and you usually want to
re-run/tune one without re-running the other.

Outputs:

  tpm_matrix.tsv    transcript_id x sample TPM, read directly from each
                    sample's GTF "TPM" attribute.
  counts_matrix.tsv the same shape, as ESTIMATED raw read counts:
                    count = round(cov * length / read_length), the standard
                    formula StringTie's own bundled prepDE.py uses to turn
                    -e/-B coverage output into a DESeq2/edgeR-compatible
                    count matrix (confirmed present at
                    <stringtie-env>/bin/prepDE.py, ships with the pinned
                    stringtie=2.2.1 bioconda package). Reimplemented inline
                    rather than shelling out to prepDE.py because that tool
                    reads a whole sample directory, is unfiltered to our
                    classified candidate set, and writes CSV rather than
                    this repo's tsv.gz convention -- this way the same
                    single pass that already extracts TPM also produces the
                    count estimate, pre-filtered and pre-formatted like
                    every other matrix here. `length` is the transcript's
                    total exonic span, summed once from --merged-gtf (the
                    same structure every sample was re-quantified against,
                    so it does not vary by sample); `read_length` is the
                    cohort's auto-detected read length (ref.smk's
                    SJDB_OVERHANG + 1) unless overridden.
  cpm_matrix.tsv    the counts matrix, each column divided by that sample's
                    column total and multiplied by 1e6 (CPM, not a second
                    TPM -- it normalizes by library size only, the same
                    convention used for the reads-screen's cpm_matrix.tsv.gz
                    in chimera_reads_counts.py). 0 for a sample whose column
                    total is 0.
"""
import argparse
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_write

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def load_candidate_ids(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("transcript_id")
        return {line.rstrip("\n").split("\t")[idx] for line in fh if line.strip()}


def load_sample_metrics(gtf_path, wanted_ids):
    """{transcript_id: {"tpm": float, "cov": float}} for transcript lines in
    `wanted_ids`."""
    metrics = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = dict(ATTR_RE.findall(cols[8]))
            tid = attrs.get("transcript_id")
            if tid in wanted_ids:
                metrics[tid] = {
                    "tpm": float(attrs.get("TPM", 0.0)),
                    "cov": float(attrs.get("cov", 0.0)),
                }
    return metrics


def load_transcript_lengths(merged_gtf_path, wanted_ids):
    """{transcript_id: total exonic span} summed from --merged-gtf's "exon"
    lines. Computed once from the shared merged structure every sample was
    re-quantified against (-e -B), not per sample -- exon boundaries for a
    given transcript_id are identical across samples by construction."""
    lengths = {}
    opener = gzip.open if merged_gtf_path.endswith(".gz") else open
    with opener(merged_gtf_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = dict(ATTR_RE.findall(cols[8]))
            tid = attrs.get("transcript_id")
            if tid not in wanted_ids:
                continue
            span = int(cols[4]) - int(cols[3]) + 1
            lengths[tid] = lengths.get(tid, 0) + span
    return lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--quant", nargs="+", required=True,
                     help="Per-sample requantified GTFs, same order as --sample-names")
    # nargs="+" like chimera_reads_counts.py / tecount_counts.py: the rule
    # interpolates the names unquoted, so the shell hands us N argv entries,
    # not one space-joined string. Without it argparse consumed only the
    # first name and exited 2 on the rest -- latent because this rule never
    # ran while chimera.assembly was off by default.
    ap.add_argument("--sample-names", required=True, nargs="+")
    ap.add_argument("--merged-gtf", required=True,
                     help="results/chimera/assembly/stringtie_merge.gtf -- "
                     "for transcript lengths (counts_matrix only)")
    ap.add_argument("--read-length", required=True, type=int,
                     help="Cohort read length (counts_matrix only)")
    ap.add_argument("--out-tpm", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-cpm", required=True)
    args = ap.parse_args()

    samples = args.sample_names
    if len(samples) != len(args.quant):
        sys.exit("--sample-names and --quant must have the same length/order")

    wanted = load_candidate_ids(args.candidates)
    lengths = load_transcript_lengths(args.merged_gtf, wanted)

    tpm_matrix = {tid: {} for tid in wanted}
    counts_matrix = {tid: {} for tid in wanted}

    for sample, gtf_path in zip(samples, args.quant):
        sample_metrics = load_sample_metrics(gtf_path, wanted)
        for tid in wanted:
            m = sample_metrics.get(tid, {"tpm": 0.0, "cov": 0.0})
            tpm_matrix[tid][sample] = m["tpm"]
            length = lengths.get(tid, 0)
            counts_matrix[tid][sample] = round(
                m["cov"] * length / args.read_length
            ) if length and args.read_length else 0

    os.makedirs(os.path.dirname(args.out_tpm) or ".", exist_ok=True)
    with open_write(args.out_tpm) as fh:
        fh.write("transcript_id\t" + "\t".join(samples) + "\n")
        for tid in sorted(tpm_matrix):
            row = [f"{tpm_matrix[tid][s]:.3f}" for s in samples]
            fh.write(tid + "\t" + "\t".join(row) + "\n")

    os.makedirs(os.path.dirname(args.out_counts) or ".", exist_ok=True)
    with open_write(args.out_counts) as fh:
        fh.write("transcript_id\t" + "\t".join(samples) + "\n")
        for tid in sorted(counts_matrix):
            row = [str(counts_matrix[tid][s]) for s in samples]
            fh.write(tid + "\t" + "\t".join(row) + "\n")

    totals = {s: sum(counts_matrix[tid][s] for tid in wanted) for s in samples}
    os.makedirs(os.path.dirname(args.out_cpm) or ".", exist_ok=True)
    with open_write(args.out_cpm) as fh:
        fh.write("transcript_id\t" + "\t".join(samples) + "\n")
        for tid in sorted(counts_matrix):
            row = [
                "0.000" if totals[s] == 0
                else f"{counts_matrix[tid][s] / totals[s] * 1e6:.3f}"
                for s in samples
            ]
            fh.write(tid + "\t" + "\t".join(row) + "\n")

    print(f"{len(wanted)} candidates x {len(samples)} samples -> "
          f"{args.out_tpm}, {args.out_counts}, {args.out_cpm}")


if __name__ == "__main__":
    main()
