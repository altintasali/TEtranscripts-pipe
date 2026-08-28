#!/usr/bin/env python3
"""Look at every line of chimera evidence side by side, instead of collapsing
it into one number.

Why this exists: chimera_evidence.tsv.gz carries several independent
measurements per (gene, TE) pair -- breakpoint reads, replicate support, a
splice motif, assembled transcripts, TE-locus expression. Each screen's
measurements are shown as its own columns, never merged into an "agreement"
flag. Reducing them to a single ordinal tier hides whether any
of them actually carries information, and twice it hid the opposite of what
was assumed: TElocal expression turned out ANTI-correlated with the splice
motif, and screen agreement turned out to sit at roughly its chance rate.
Both were obvious the moment the dimensions were plotted against each other
and invisible while they were being summed.

So this emits two heatmaps and no score:

  correlation   dimension x dimension (Spearman) across all pairs the
                junction screen found. Reads as: which of these say the same
                thing (redundant), which are independent, and which point
                the wrong way. A dimension that correlates with nothing is
                either noise or the only real information here -- the rest
                of the report has to decide which.

  candidates    the union of the top-N pairs on EACH dimension, scored as a
                within-cohort percentile so the columns are comparable.
                Deliberately not "the top N by some combined score": the
                point is to show whether the leaders on one axis lead on any
                other. A blocky, mostly-diagonal picture means they do not,
                and that no single ranking is defensible yet.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read

# (column, short label, how to read a row's value)
#
# Each screen contributes its own columns and there is deliberately NO
# combined "both screens" column. It used to be here, derived as
# (junction_events > 0 and assembly_transcripts > 0) -- but it is a function
# of the columns either side of it, so it added no information while looking
# like an independent line of evidence. Collapsing two screens into one
# boolean is the same premature aggregation that hid the tier-4 problem: read
# the two screens' columns and judge the agreement yourself.
DIMENSIONS = [
    ("junction_reads", "Junction reads", "num"),
    ("junction_events", "Breakpoints", "num"),
    ("junction_max_samples", "Replicates", "num"),
    ("junction_canonical", "Splice motif", "yes"),
    ("assembly_transcripts", "Assembly transcripts", "num"),
    ("assembly_strand_match", "Assembly strand match", "yes"),
    ("telocal_active", "TE locus expressed", "yes"),
]


def value(row, col, kind):
    if kind == "num":
        try:
            return float(row.get(col, 0) or 0)
        except ValueError:
            return 0.0
    if kind == "yes":
        return 1.0 if row.get(col) == "yes" else 0.0
    return 0.0


def ranks(values):
    """Average ranks, ties shared -- the basis for Spearman without scipy."""
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = shared
        i = j + 1
    return out


def pearson(a, b):
    n = len(a)
    if n < 2:
        return 0.0
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    if va <= 0 or vb <= 0:
        # A dimension with no variation (e.g. every pair canonical, or the
        # assembly screen disabled) has no correlation to report. 0 is the
        # honest value; it must not be mistaken for "measured, unrelated".
        return 0.0
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / (va ** 0.5 * vb ** 0.5)


def load(path):
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def heatmap_doc(doc_id, section, description, xcats, ycats, data, title):
    return {
        "id": doc_id,
        # Same group as every other view of the read-evidence screen.
        "parent_id": "chimera_reads",
        "parent_name": "Gene-TE chimeras: read evidence",
        "section_name": section,
        "description": description,
        "plot_type": "heatmap",
        "pconfig": {"id": f"{doc_id}_plot", "title": title},
        "xcats": xcats,
        "ycats": ycats,
        "data": data,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--gene-names", default=None,
                    help="gene_id_to_name.tsv.gz, to label rows by symbol")
    ap.add_argument("--out-correlation", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--top-n", type=int, default=8,
                    help="per dimension; the candidate heatmap is their union")
    args = ap.parse_args()

    rows = list(load(args.evidence))

    symbols = {}
    if args.gene_names and os.path.exists(args.gene_names):
        with open_read(args.gene_names) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    symbols[f[0]] = f[1]

    labels = [lbl for _, lbl, _ in DIMENSIONS]

    # --- correlation ----------------------------------------------------
    # Restricted to pairs the junction screen found. Assembly-only pairs have
    # every junction_* field at 0, which would manufacture correlation between
    # all of them purely from that shared absence.
    jrows = [r for r in rows if r.get("found_by") in ("junction", "both")]
    cols = [[value(r, c, k) for r in jrows] for c, _, k in DIMENSIONS]
    rcols = [ranks(c) for c in cols]
    matrix = [[round(pearson(rcols[i], rcols[j]), 3) for j in range(len(cols))]
              for i in range(len(cols))]

    corr_desc = (
        f"Spearman correlation between the evidence columns of "
        f"chimera_evidence.tsv.gz, across the {len(jrows):,} pair(s) the "
        "junction screen found (assembly-only pairs are excluded: their "
        "junction fields are all zero, which would manufacture correlation "
        "from shared absence). "
        "<br><br><em>How to read it:</em> two dimensions that correlate "
        "strongly are telling you the same thing, so treating them as "
        "separate support double-counts. A dimension that correlates with "
        "nothing else is either noise or the only independent information in "
        "the table. A <strong>negative</strong> correlation with "
        "<em>Splice motif</em> is the one to take seriously: the motif is "
        "the best artifact discriminator here, so anything anti-correlated "
        "with it is selecting against real splicing, not for it."
    )
    corr = heatmap_doc(
        "chimera_evidence_correlation", "Evidence correlation", corr_desc,
        labels, labels, matrix, "Spearman correlation between evidence types",
    )

    # --- candidates -----------------------------------------------------
    # Percentile within the cohort, so columns on wildly different scales
    # (1-2 reads vs 14 transcripts) can sit side by side.
    pct = []
    n = len(rows)
    for c, _, k in DIMENSIONS:
        vals = [value(r, c, k) for r in rows]
        rk = ranks(vals)
        pct.append([round(100.0 * x / n, 1) for x in rk] if n else [])

    picked, seen = [], set()
    for di in range(len(DIMENSIONS)):
        order = sorted(range(n), key=lambda i: pct[di][i], reverse=True)
        for i in order[: args.top_n]:
            if i not in seen:
                seen.add(i)
                picked.append(i)

    def label(i):
        r = rows[i]
        gene = symbols.get(r.get("gene_id", ""), r.get("gene_id", "?"))
        return f"{gene} / {r.get('te_id', '?')}"

    cand_desc = (
        f"The union of the top {args.top_n} pair(s) on <em>each</em> evidence "
        "dimension, so a row is here because it leads on at least one axis. "
        "Cells are within-cohort percentiles, not raw values, so columns "
        "measured on different scales are comparable. "
        "<br><br><em>How to read it:</em> a row that is bright across many "
        "columns has broad support. A row bright in one column and dark in "
        "the rest leads on that axis alone -- which is most of them, on real "
        "data. If the picture is mostly diagonal, the dimensions disagree "
        "about what the best candidates are, and no single ranking of this "
        "table is defensible yet. That is a finding, not a failure of the "
        "plot: it says the evidence is genuinely thin and the candidates "
        "should be judged individually."
    )
    cand = heatmap_doc(
        "chimera_evidence_candidates", "Candidates by evidence type",
        cand_desc, labels, [label(i) for i in picked],
        [[pct[d][i] for d in range(len(DIMENSIONS))] for i in picked],
        "Evidence percentile per candidate",
    )

    for path, doc in ((args.out_correlation, corr), (args.out_candidates, cand)):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")

    print(f"evidence heatmaps: {len(rows)} pairs ({len(jrows)} junction-side), "
          f"{len(picked)} candidates shown")


if __name__ == "__main__":
    main()
