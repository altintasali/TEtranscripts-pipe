#!/usr/bin/env python3
"""Per-sample TEcounts summary stats for the MultiQC report (custom content):

TEcount writes one {sample}.cntTable per sample: a two-column table
(`gene/TE` feature key + read count) holding gene ids and TE subfamily names
mixed in one column. TEcount builds TE keys as `gene_id:family_id:class_id`
from the TE GTF, so each row can be classified from the key alone: a key with
no colons is a gene; a `gene:family:class` key is a TE subfamily whose family
and class come from the key.

Emits:
  tecount_assignment_mqc.json  genes vs TE subfamilies (counts + % of total)
  tecount_te_class_mqc.json    TE subfamily reads by class, LINE/SINE/LTR/
                               DNA/RC + unknown

Unlike the sample-QC matrix (tecount_counts.py), this uses the RAW cntTables
(all features), so it is independent of tetranscripts.qc.feature_class.

The implementation is shared with telocal_summary_mqc.py (te_summary_common):
the two sections answer the same question at two resolutions, and keeping two
copies let them drift -- see that module's docstring for the silent-zero bug
that came of it.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from te_summary_common import Flavour, run

FLAVOUR = Flavour(
    tool="TEcount",
    # Must match sample_qc.R's tecount view id -- a mismatch splits the
    # TEcount group into two report sections.
    parent_id="tecount",
    unit_label="TE subfamilies",
    # gene_id:family_id:class_id -- 3 fields, but 2 is still a TE key.
    min_key_fields=2,
    assignment_title="TEcounts assignment (genes vs TEs)",
    class_title="TE class composition",
    assignment_desc=(
        "Per-sample read counts assigned by TEcount to genes vs TE "
        "subfamilies (counts and % of total assigned reads). "
        "<br><br><em>How to read this:</em> the gene/TE split should be "
        "broadly consistent across samples of the same type — an "
        "outlier usually reflects library quality or rRNA/intronic "
        "carry-over rather than TE biology, so check it before "
        "interpreting a TE effect. TEcount pools every copy of a "
        "subfamily into one row; see the TElocal sections for the "
        "per-copy view."
    ),
    class_desc=(
        "Per-sample TE-subfamily reads by repeat class (counts and % of "
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
