"""Turn a per-sample splice-junction table (classify_chimera_splice_junctions.py)
into a BED track for IGV: one BED6-ish row per gene-TE junction, spanning the
intron itself. Mirrors chimera_chimeric_reads_to_igv_bed.py's shape exactly,
so the two screens' tracks load and read the same way side by side in IGV.

The script runs under Snakemake's `script:` directive, so it reads the
snakemake.input / snakemake.output / snakemake.params globals.

BED columns written (BED6, score = unique reads, strand = junction strand):
    chrom  intron_start  intron_end  event_id  unique_reads  strand
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read

with open_read(snakemake.input[0]) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    rows = [
        dict(zip(header, line.rstrip("\n").split("\t")))
        for line in fh
        if line.strip()
    ]

os.makedirs(os.path.dirname(str(snakemake.output[0])), exist_ok=True)
with open(snakemake.output[0], "w") as fh:
    fh.write('track name="chimera_splice_junctions" description="gene-TE '
             'splice junctions ({sample})" itemRgb="On"\n')
    for r in rows:
        if r.get("direction") not in ("gene_to_te", "te_to_gene"):
            continue
        try:
            start, end = int(r["intron_start"]), int(r["intron_end"])
        except (KeyError, ValueError):
            continue
        # BED is 0-based half-open; SJ.out.tab intron coordinates are
        # 1-based inclusive (first/last base of the intron), so the intron
        # itself spans [start-1, end).
        fh.write("\t".join([
            r.get("chrom", "."), str(start - 1), str(end),
            r.get("event_id", "."), r.get("unique_reads", "1"),
            r.get("strand", "."),
        ]) + "\n")
