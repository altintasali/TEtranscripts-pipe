#!/usr/bin/env python3
"""The report's gene-TE candidate list, as a table the reader sorts.

This is the entry point users actually want: a look at the candidates without
opening a TSV. What it deliberately is NOT is a ranking.

The pipeline has twice had a ranking removed. A four-tier confidence ladder
went first, for having no validated weighting; three competing per-screen
top-N tables went with it. Re-introducing "ordered by splice motif, then
replicates, then strand match, then depth" would be that same ordering under a
disclosure line -- it is, key for key, the rank_key deleted in d927c8f.

So the ordering is handed to the reader instead. MultiQC's native table
(plot_type: "table") sorts on any column, every signal is its own column, and
the default order is n_evidence -- a COUNT of how many evidence types a pair
carries, which weights nothing precisely because no weighting has been
established. Which column deserves weight is what the guide section above
explains; the reader applies it by clicking a header.

Row cap: a real cohort produces tens of thousands of pairs and MultiQC embeds
table data in the HTML, so only the top --top-n are rendered. The full
catalogue is candidates.tsv.gz, and the section says so.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

PARENT_ID = "chimera"
PARENT_NAME = "Chimera"

# (column in candidates.tsv.gz, header, description, kind)
# kind: "int" -> numeric, right-aligned; "yesno" -> rendered as a yes/no string
COLUMNS = [
    # Gene and TE get their own columns, not just the row key: the key is one
    # string, so without these you cannot sort by gene, and a reader scanning
    # for a specific gene has nothing to sort on.
    ("gene_label", "Gene",
     "Gene symbol where the reference GTF provides one, otherwise the "
     "gene_id.", "str"),
    ("te_id", "TE insertion",
     "The individual TE copy (transcript_id in the TE GTF), not the "
     "subfamily. Joins against TElocal rows.", "str"),
    ("n_evidence", "Evidence types",
     "How many of the four evidence flags this pair carries. A count, not a "
     "score -- the flags are unweighted.", "int"),
    ("junction_canonical", "Splice motif",
     "A recognised splice motif on at least one junction (STAR). The guide "
     "above calls this the best artifact discriminator available.", "yesno"),
    ("junction_max_samples", "Samples",
     "Most samples any one junction for this pair was seen in (STAR).", "int"),
    ("found_by", "Found by",
     "Which screens called it: junction (STAR), assembly (StringTie), or "
     "both. Agreement measured near its chance rate -- see the guide.", "str"),
    ("assembly_strand_match", "Strand match",
     "The assembled transcript's strand agrees with the gene's (StringTie).",
     "yesno"),
    ("junction_reads", "Reads",
     "Chimeric reads supporting this pair (STAR). The metric most inflated by "
     "artifacts -- shown last on purpose.", "int"),
    # Reported, never counted as evidence. It is in candidates.tsv.gz and the
    # guide discusses it at length, so leaving it out of the table meant the
    # one place a reader looks did not show it -- and its absence read as the
    # column not existing rather than as a deliberate exclusion.
    ("telocal_active", "TE locus expressed",
     "TElocal: whether the TE copy is itself transcribed. CONTEXT, NOT "
     "EVIDENCE -- it does not contribute to Evidence types. Measured on a "
     "real cohort it was anti-correlated with the splice motif (6.7% vs "
     "10.2% canonical), so an expressed locus is not support. '.' means "
     "TElocal did not run.", "str"),
]


def load(path):
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def load_symbols(path):
    symbols = {}
    if not path or not os.path.exists(path):
        return symbols
    with open_read(path) as fh:
        fh.readline()
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1] not in (".", ""):
                symbols[fields[0]] = fields[1]
    return symbols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True,
                    help="results/chimera/candidates.tsv.gz")
    ap.add_argument("--gene-names", default=None,
                    help="gene_id_to_name.tsv.gz, to label rows by symbol")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--source-path", default="results/chimera/candidates.tsv.gz")
    args = ap.parse_args()

    rows = list(load(args.evidence))
    symbols = load_symbols(args.gene_names)

    # candidates.tsv.gz is already written in n_evidence order; take the head
    # rather than re-sorting, so the report and the file agree exactly.
    top = rows[: args.top_n]

    data = {}
    for r in top:
        gene_id = r.get("gene_id", ".")
        gene = symbols.get(gene_id, gene_id)
        # " | ", never " / ": MultiQC cleans table row names like filenames
        # and splits on "/", taking the basename -- which silently dropped the
        # gene and left only the TE id. Measured; guard 50 pins it.
        key = f"{gene} | {r.get('te_id', '.')}"
        # a duplicate key would silently drop a row; disambiguate with gene_id
        if key in data:
            key = f"{key} ({gene_id})"
        entry = {}
        for col, header, _desc, kind in COLUMNS:
            value = gene if col == "gene_label" else r.get(col, ".")
            entry[header] = _int(value) if kind == "int" else str(value)
        data[key] = entry

    headers = {
        header: {
            "title": header,
            "description": desc,
            **({"format": "{:,.0f}", "min": 0} if kind == "int" else {}),
        }
        for _col, header, desc, kind in COLUMNS
    }

    note = (
        "<p>The <strong>{n_shown:,}</strong> of {n_total:,} gene-TE pairs with "
        "the most evidence types, from <code>{src}</code>. "
        "<strong>Click any column header to sort.</strong></p>"
        "<p>The default order is <em>Evidence types</em>, a count of how many "
        "of the four flags a pair carries. It is not a score: the flags are "
        "counted unweighted, because no weighting of them has been validated "
        "here. It also tilts toward pairs the assembly screen found, since two "
        "of the four flags need assembly support. <strong>How to weigh this "
        "evidence</strong> above says what each column is actually worth; "
        "sort on the one your question needs, and validate candidates "
        "manually before treating any of them as a result.</p>"
    ).format(n_shown=len(top), n_total=len(rows), src=args.source_path)

    if top:
        body = {
            "plot_type": "table",
            "pconfig": {
                "id": "chimera_candidates_table_plot",
                "title": "Gene-TE chimera candidates",
                "col1_header": "Gene | TE insertion",
                # defaultsort is what actually sets the opening order.
                # sort_rows: False does NOT survive -- MultiQC re-populates its
                # camelCase alias sortRows from the default (True), so the rows
                # arrive alphabetised by name whatever this says. Measured.
                # Stating the intended sort explicitly is the reliable route.
                "sort_rows": False,
                "defaultsort": [{"column": "Evidence types", "direction": "desc"}],
                "no_violin": True,
            },
            "headers": headers,
            "data": data,
        }
    else:
        # An empty table is a normal outcome; MultiQC crashes the whole report
        # on a plot with no data, so fall back to prose (see
        # chimera_evidence_guide_mqc.py for the same trap).
        body = {
            "plot_type": "html",
            "data": (
                "<p>No gene-TE chimera candidates were found in this run. "
                "This is not an error: it means neither screen made a call. "
                "TEcount and TElocal results are unaffected.</p>"
            ),
        }

    doc = {
        # must NOT be "chimera": a section id equal to a parent_id is picked up
        # by report_section_order's module pass and sends the group to the end.
        "id": "chimera_candidates_table",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "section_name": "Candidates",
        "description": note,
        **body,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"chimera candidates table: {len(top)} of {len(rows)} pairs shown "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
