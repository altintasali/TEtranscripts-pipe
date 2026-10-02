#!/usr/bin/env python3
"""Annotated splice features of the gene GTF, for telling KNOWN gene
structure apart from new gene-TE chimeras.

Two outputs:

  annotated_introns.tsv.gz  every intron of every annotated transcript, as
                            chrom, intron_start, intron_end, strand --
                            1-based and inclusive, i.e. exactly STAR's
                            SJ.out.tab columns 1-3, so an SJ junction is
                            annotated iff (chrom, start, end) is in this
                            set. Used by classify_chimera_splice_junctions.py
                            to type such junctions annotated_splice.
                            STAR's own "annotated" column is NOT used for
                            this: under star.two_pass: cohort STAR marks
                            every junction it inserted from the pooled
                            first pass as annotated too (on a real run 51%
                            of chimera junctions by STAR's flag vs 25%
                            against the GTF).
  last_exons.bed            BED6, one row per annotated transcript's last
                            exon (3'-most, strand-aware): chrom, start,
                            end, transcript_id, score, strand -- the 3'
                            counterpart of annotation_to_bed.py's
                            first_exons.bed. Used by
                            classify_chimera_assembly.py to type a
                            TE-overlapping last exon that is an annotated
                            last exon (a TE in an ordinary 3' UTR)
                            annotated_terminal_exon_embedded_te instead of
                            te_terminated (43% of te_terminated calls on a
                            real run).

Kept separate from annotation_to_bed.py on purpose: changing that rule
would rewrite genes.bed/exons.bed/te.bed and re-run everything built on
them. This rule only feeds the two classifiers above.

Usage:
    annotation_splice_features.py --gtf GENE.GTF \\
        --out-introns annotated_introns.tsv.gz --out-last-exons last_exons.bed
"""
import argparse
import os
import sys
from itertools import pairwise

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotation_to_bed import parse_attrs
from gz_io import open_read, open_write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-introns", required=True)
    ap.add_argument("--out-last-exons", required=True)
    args = ap.parse_args()

    transcripts = {}
    with open_read(args.gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            try:
                s, e = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            a = parse_attrs(cols[8])
            tid = a.get("transcript_id") or a.get("gene_id")
            if not tid:
                continue
            t = transcripts.setdefault(
                tid, {"chrom": cols[0], "strand": cols[6], "exons": []})
            t["exons"].append((s, e))

    introns = set()
    last_rows = []
    for tid, t in transcripts.items():
        ex = sorted(t["exons"])
        for (_s1, e1), (s2, _e2) in pairwise(ex):
            if s2 - 1 >= e1 + 1:
                introns.add((t["chrom"], e1 + 1, s2 - 1, t["strand"]))
        last_s, last_e = ex[0] if t["strand"] == "-" else ex[-1]
        last_rows.append((t["chrom"], last_s - 1, last_e, tid, ".", t["strand"]))

    for path in (args.out_introns, args.out_last_exons):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open_write(args.out_introns) as fh:
        fh.write("chrom\tintron_start\tintron_end\tstrand\n")
        for row in sorted(introns):
            fh.write("\t".join(str(x) for x in row) + "\n")
    with open_write(args.out_last_exons) as fh:
        for row in sorted(last_rows):
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"wrote {len(introns)} annotated introns and {len(last_rows)} "
          f"last exons from {len(transcripts)} transcripts")


if __name__ == "__main__":
    main()
