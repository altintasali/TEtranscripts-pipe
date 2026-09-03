#!/usr/bin/env python3
"""Build a transcript_id x sample TPM matrix for chimera-assembly candidates,
pulling values out of each sample's -e/-B re-quantified StringTie GTF
(results/chimera/assembly/per_sample/quant/{sample}.transcripts.gtf, see
chimera_assembly.smk's stringtie_requantify rule).

Kept separate from classify_chimera_assembly.py the same way chimera_reads_counts.py
is separate from classify_chimera_reads.py: structural classification and
expression aggregation are different concerns, and you usually want to
re-run/tune one without re-running the other.
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


def load_sample_tpm(gtf_path, wanted_ids):
    """{transcript_id: tpm} for transcript lines in `wanted_ids`."""
    tpm = {}
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
                tpm[tid] = float(attrs.get("TPM", 0.0))
    return tpm


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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    samples = args.sample_names
    if len(samples) != len(args.quant):
        sys.exit("--sample-names and --quant must have the same length/order")

    wanted = load_candidate_ids(args.candidates)
    matrix = {tid: {} for tid in wanted}

    for sample, gtf_path in zip(samples, args.quant):
        sample_tpm = load_sample_tpm(gtf_path, wanted)
        for tid in wanted:
            matrix[tid][sample] = sample_tpm.get(tid, 0.0)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("transcript_id\t" + "\t".join(samples) + "\n")
        for tid in sorted(matrix):
            row = [f"{matrix[tid][s]:.3f}" for s in samples]
            fh.write(tid + "\t" + "\t".join(row) + "\n")

    print(f"TPM matrix: {len(matrix)} candidates x {len(samples)} samples -> {args.out}")


if __name__ == "__main__":
    main()
