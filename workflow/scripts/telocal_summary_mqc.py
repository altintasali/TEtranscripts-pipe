#!/usr/bin/env python3
"""Per-sample TElocal summary stats for the MultiQC report (custom content):

TElocal writes one {sample}.cntTable per sample: a two-column table
(`gene/TE` feature key + read count) holding gene ids and TE locus keys
mixed in one column. TElocal builds TE keys as
`transcript_id:gene_id:family_id:class_id` from the .locInd index, where
transcript_id embeds genomic coordinates (e.g.
`chr1:564318:564741(L1PA2:+):L1PA2:L1:LINE/L1`). Each row can be classified
from the key alone: a key with no colons is a gene; a key with >=2 colons
has family = second-to-last field, class = last field.

Emits:
  telocal_assignment_mqc.json  genes vs TE loci (counts + % of total)
  telocal_te_class_mqc.json    TE locus reads by class, LINE/SINE/LTR/
                               DNA/RC + unknown

Unlike the sample-QC matrix, this uses the RAW cntTables (all features).

The implementation is shared with tecount_summary_mqc.py (te_summary_common):
the two sections answer the same question at two resolutions, and keeping two
copies let them drift -- this file's own count parser used to turn a
float-formatted count into a silent zero.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from te_summary_common import Flavour, run

FLAVOUR = Flavour(
    tool="TElocal",
    # Must match sample_qc.R's telocal view id.
    parent_id="telocal",
    unit_label="TE loci",
    # transcript_id:gene_id:family_id:class_id -- and transcript_id itself
    # contains colons, so a TE key never has fewer than 3 fields.
    min_key_fields=3,
    assignment_title="TElocal assignment (genes vs TEs)",
    class_title="TElocal TE class composition",
    assignment_desc=(
        "Per-sample read counts assigned by TElocal to genes vs TE "
        "loci (counts and % of total assigned reads). "
        "<br><br><em>How to read this:</em> TElocal resolves each TE "
        "insertion separately, where TEcount pools all copies of a "
        "subfamily into a single row. A subfamily that looks uniformly "
        "expressed in TEcount is often one or two highly-expressed "
        "copies here, with the rest silent — which is the reason "
        "to run both."
    ),
    class_desc=(
        "Per-sample TE-locus reads by repeat class (counts and % of "
        "TE reads); classes not in LINE/SINE/LTR/DNA/RC are grouped as "
        "unknown."
    ),
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--out-assignment", required=True)
    ap.add_argument("--out-class", required=True)
    run(FLAVOUR, ap.parse_args())


if __name__ == "__main__":
    main()
