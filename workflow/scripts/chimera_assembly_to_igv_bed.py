#!/usr/bin/env python3
"""Turn results/chimera_assembly/candidates.tsv.gz (classify_chimera_assembly.py)
into a BED track for IGV: one row per candidate, spanning the SPECIFIC exon
that overlaps the TE (te_exon_start/te_exon_end -- not the whole transcript,
and not always the first exon; see classify_chimera_assembly.py's docstring),
colored by chimera_type so the different classes are visually distinguishable
when loaded alongside the chimera_junction screen's own IGV track.

candidates.tsv.gz coordinates are already GTF-derived 0-based-start/
exclusive-end (see classify_chimera_assembly.py's load_transcripts), i.e.
already valid BED coordinates -- no conversion needed here.

BED columns written (BED9): chrom, te_exon_start, te_exon_end, name, score,
strand, thickStart, thickEnd, itemRgb.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read

CLASS_COLOR = {
    "te_initiated": "31,119,180",
    "te_initiated_intergenic": "150,190,220",
    "te_exonized": "255,127,14",
    "te_terminated": "44,160,44",
    "unspliced_te_only": "150,150,150",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True,
                     help="results/chimera_assembly/candidates.tsv.gz (or "
                          "the _with_junction_evidence variant)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open_read(args.candidates) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n = 0
    with open(args.out, "w") as fh:
        fh.write('track name="chimera_assembly" description="StringTie-assembly '
                  'chimera candidates" itemRgb="On"\n')
        for r in rows:
            try:
                start, end = int(r["te_exon_start"]), int(r["te_exon_end"])
            except (KeyError, ValueError):
                continue
            chrom = r.get("chrom", ".")
            strand = r.get("strand", ".")
            chimera_type = r.get("chimera_type", ".")
            name = f"{r.get('transcript_id', '.')};{chimera_type}"
            color = CLASS_COLOR.get(chimera_type, "0,0,0")
            score = "900" if r.get("confirmed_by_junction_screen") == "yes" else "300"
            fh.write("\t".join([
                chrom, str(start), str(end), name, score, strand,
                str(start), str(end), color,
            ]) + "\n")
            n += 1

    print(f"{n} chimera-assembly candidates -> {args.out}")


if __name__ == "__main__":
    main()
