#!/usr/bin/env python3
"""Per-locus genomic coordinates for TElocal's TE keys, as BED6.

TElocal's cntTable reports TE counts keyed by
`transcript_id:gene_id:family_id:class_id` (see build_telocal_index.py and
telocal_summary_mqc.py), but the cntTable itself never includes genomic
coordinates -- that's how TElocal's own output is shaped. This script does
one pass over the same TE GTF used to build the .locInd index and emits a
standard BED6 file (chrom, start, end, name, score, strand) with that key as
the `name` column, so cntTable rows can be joined against real genomic
coordinates (e.g. to plot or intersect with other data) while also getting
native bedtools/IGV compatibility for free -- joining on a named column
costs nothing over joining on column 1. No index/tree needed -- it just
reads the GTF's own coordinate columns, so it works regardless of what
naming convention a given TE GTF happens to use for transcript_id.

Follows the same BED conventions as annotation_to_bed.py's genes.bed/
exons.bed/te.bed: 0-based half-open start, "." for the unused score field,
no header line, uncompressed (plain gzip isn't random-access-indexed the
way IGV needs -- it wants plain text or bgzip+tabix -- so gzipping here
would break the exact IGV-loadable property this format switch is for).

Usage:
  telocal_locations.py --gtf te.gtf --out results/telocal/telocal_locations.bed
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_telocal_index import _ATTR_RE
from gz_io import open_read, open_write


def build_locations(gtf_path, out_path):
    seen = set()
    with open_read(gtf_path) as fh, open_write(out_path) as out:
        linenum = 0
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split("\t")
            chrom = items[0]
            start0 = int(items[3]) - 1  # GTF is 1-based inclusive; BED start is 0-based
            end = items[4]
            strand = items[6]
            linenum += 1

            attrs = dict(_ATTR_RE.findall(items[8]))
            ele_id = attrs.get("gene_id", "")
            transc_id = attrs.get("transcript_id", "")
            family_id = attrs.get("family_id", "")
            class_id = attrs.get("class_id", "")
            if ele_id == "" or transc_id == "" or family_id == "" or class_id == "":
                sys.stderr.write(line + "\n")
                sys.stderr.write(
                    f"TE GTF format error! There is no annotation at line {linenum}.\n"
                )
                raise ValueError(f"TE GTF format error at line {linenum}")

            # Matches build_telocal_index.py's ele_name exactly -- the same
            # string TElocal's cntTable uses as its "gene/TE" row key.
            ele_name = f"{transc_id}:{ele_id}:{family_id}:{class_id}"
            if ele_name in seen:
                # First occurrence wins, same as _elements/.locInd -- TE
                # GTFs are expected to give every instance a unique key.
                continue
            seen.add(ele_name)
            out.write(f"{chrom}\t{start0}\t{end}\t{ele_name}\t.\t{strand}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Emit a BED6 of chrom/start/end/strand for every TElocal TE key."
    )
    ap.add_argument("--gtf", required=True, help="TE GTF (plain or gzipped)")
    ap.add_argument(
        "--out", required=True, help="Output BED6 path (plain text, for IGV/bedtools)"
    )
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    build_locations(args.gtf, args.out)


if __name__ == "__main__":
    main()
