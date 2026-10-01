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
the default order is Screens descending only -- a COUNT of how many
independent screens found a pair, which weighs nothing precisely because no
weighting has been established. Two per-pair counts used to break ties here
before (n_screen_evidence, n_corroboration); they are gone from this view
because the flags behind them are now their own columns (CR/SJ motif,
Assembly strand, Replicated, TElocal reads) -- summarising them a second
time as a silent tie-break was exactly the redundancy this rewrite removes.
Which column deserves weight is what the guide section below explains; the
reader applies it by clicking a header.

Columns are grouped by screen and coloured accordingly (CR = blue, SJ =
green, Assembly = orange, Support = purple) so the blocks are visible at a
glance; Gene and TE insertion are still their own columns (so the reader can
sort/filter on either alone) but hidden by default since both already
appear in the row key.

Row selection: a real cohort produces tens of thousands of pairs and MultiQC
embeds table data in the HTML, so this shows every pair tied at the top
Screens value (filling down to MIN_ROWS from the next Screens value when
that group is small), capped at --max-rows as a hard safety ceiling -- see
_select_rows()'s docstring. The full catalogue is candidates.tsv.gz, and the
section says so.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write
from chimera_exon_context import CHIMERA_CALL_TYPES

PARENT_ID = "chimera"
PARENT_NAME = "Chimera"

# (column in candidates.tsv.gz, header, description, kind)
# kind: "int" -> numeric, right-aligned; "yesno" -> rendered as a yes/no
# string; "replicated" is the one synthetic column -- not a candidates.tsv.gz
# field, derived in the loop below from the corroboration column instead.
#
# Grouped by block, and within a block always motif -> samples -> reads, so
# the three screens read the same way left to right: Pair, TE vs gene,
# Overview (cross-screen), then one block per screen (CR / SJ / Assembly), then
# Support (cross-cutting corroboration). The CR / SJ / Assembly header
# prefixes deliberately match found_by's own tokens (cr / sj / assembly) --
# see that column's description.
COLUMNS = [
    # Gene and TE get their own columns, not just the row key: the key is one
    # string, so without these you cannot sort by gene, and a reader scanning
    # for a specific gene has nothing to sort on. Hidden by default (see
    # HIDDEN_HEADERS below) since both already appear in the row key -- kept
    # as real columns, not dropped, so sort/filter-by-gene still works.
    ("gene_label", "Gene",
     "Gene symbol where the reference GTF provides one, otherwise the "
     "gene_id. Hidden by default -- already shown in the row key; use this "
     "column to sort or filter on gene alone.", "str"),
    ("te_id", "TE insertion",
     "The individual TE copy (transcript_id in the TE GTF), not the "
     "subfamily. Joins against TElocal rows. Hidden by default -- already "
     "shown in the row key.", "str"),
    ("te_subfamily", "TE subfamily",
     "TE annotation field from the curated TE GTF.", "str"),
    ("te_class", "TE class",
     "TE annotation field from the curated TE GTF.", "str"),
    # TE vs gene: where the TE sits relative to its gene and on which
    # strand, from the annotation alone (chimera_evidence.py).
    ("te_position", "TE position",
     "Where the TE sits relative to the gene, strand-aware: upstream (5' "
     "of the gene), intronic, exonic or downstream (3'). From the "
     "annotation, not from the reads.", "str"),
    ("te_orientation", "TE orientation",
     "sense / antisense: the TE's annotated strand vs the gene's. Needs no "
     "stranded library. Read it only together with TE position -- each "
     "position has its own background mix; see the guide.", "str"),
    ("te_gene_distance_bp", "Distance",
     "Gap in bp between the TE and the gene's annotated span; 0 when the "
     "TE overlaps the gene (intronic / exonic).", "intna"),
    ("n_screens", "Screens",
     "How many of the 3 independent detection screens (cr, assembly, sj) "
     "found this pair (1-3). Informational, derived from Found by.", "int"),
    ("found_by", "Found by",
     "Which screens called it, \"+\"-joined: cr (chimeric-reads screen, "
     "STAR), assembly (StringTie), sj (SJ.out.tab screen, STAR) -- e.g. "
     "\"cr+sj\". How much weight agreement between screens deserves has "
     "not been established -- see the guide. The CR / SJ / Assembly "
     "column-block headers below match these same three tokens.", "str"),
    # Synthetic, like replicated below: the union of the chimera calls any
    # screen made for the pair, from the three *_chimera_types columns.
    ("types", "Types",
     "Chimera calls made for this pair by any screen, comma-joined: "
     "te_initiated, te_terminated, te_exonized. Other types a screen "
     "reported (antisense_to_gene, known gene structure) are not calls and "
     "are listed per screen in candidates.tsv.gz.", "str"),
    ("cr_canonical", "CR motif",
     "A recognised splice motif on at least one chimeric-junction read "
     "(STAR, chimeric-reads screen). The guide below calls this the best "
     "artifact discriminator available.", "yesno"),
    ("cr_max_samples", "CR samples",
     "Most samples any one chimeric junction for this pair was seen in "
     "(STAR, chimeric-reads screen).", "int"),
    ("cr_reads", "CR reads",
     "Chimeric reads supporting this pair (STAR, chimeric-reads screen), "
     "summed across every sample that saw any of this pair's junction "
     "events -- a real cohort total. The metric most inflated by "
     "artifacts -- shown last in this block on purpose.", "int"),
    ("sj_canonical", "SJ motif",
     "A recognised splice motif on at least one SJ.out.tab junction "
     "(chimera.splice_junctions screen) -- kept separate from CR motif "
     "above: a structurally independent measurement (STAR's normal-splice "
     "motif call, not the chimeric-junction one). \".\" when that screen "
     "is disabled.", "yesno"),
    ("sj_max_samples", "SJ samples",
     "Most samples any one SJ.out.tab junction for this pair was seen in "
     "(chimera.splice_junctions screen).", "int"),
    ("sj_reads", "SJ reads",
     "STAR-reported unique reads across this pair's SJ.out.tab junctions "
     "(chimera.splice_junctions screen), summed across every sample that "
     "saw any of them -- a real cohort total, same caveat as CR reads "
     "above.", "int"),
    ("assembly_strand_match", "Assembly strand",
     "The assembled transcript's strand agrees with the gene's (StringTie, "
     "assembly screen).", "yesno"),
    ("assembly_transcripts", "Assembly transcripts",
     "Number of StringTie-assembled transcripts classified as this "
     "gene-TE chimera (StringTie, assembly screen).", "int"),
    # Cross-cutting corroboration -- can fire regardless of how many screens
    # found the pair, unlike the screen-bound flags above. Not a
    # candidates.tsv.gz column: derived here from the corroboration column,
    # mirroring chimera_evidence.py's own multi_sample rule exactly (see the
    # data-population loop below) rather than re-deriving it independently.
    ("replicated", "Replicated",
     "Seen in more than one sample by either the CR or SJ screen (the "
     "greater of CR samples / SJ samples is above 1). Corroboration, not "
     "screen-bound evidence -- can fire from a single screen alone. "
     "Mirrors the multi_sample flag in candidates.tsv.gz's corroboration "
     "column.", "yesno"),
    # Reported, never counted as evidence. It is in candidates.tsv.gz and the
    # guide discusses it at length, so leaving it out of the table meant the
    # one place a reader looks did not show it -- and its absence read as the
    # column not existing rather than as a deliberate exclusion.
    ("telocal_count", "TElocal reads",
     "TElocal read count for the TE copy itself, summed across every sample "
     "(chimera_telocal_index.py builds one shared, all-samples index up "
     "front, so this is already a cohort total, not a single sample's "
     "count -- candidates_explorer.html's own \"TElocal reads\" column "
     "reaches the same number by joining "
     "counts_matrix.tsv.gz directly). Counted as corroboration when "
     "nonzero -- standing is not validated, see the guide. Blank means "
     "TElocal did not run; 0 means it ran and found nothing.", "intna"),
]

# Columns present so the reader can sort/filter on gene or TE alone, but
# hidden from the initial view since both already appear in the row key --
# same "still there, just not always-rendered" pattern used for the guide's
# detailed table (moved into helptext, not deleted).
HIDDEN_HEADERS = {"Gene", "TE insertion"}

# Colour each screen's block distinctly (only the numeric columns -- a
# background colour scale reduces a yes/no string to a single flat colour
# via MultiQC's own scale lookup, which is not useful, so those columns keep
# the table's plain default instead). Verified by rendering: MultiQC 1.33
# accepts a named ColorBrewer scale per header with no other config.
SCALE_BY_HEADER = {
    "CR samples": "Blues", "CR reads": "Blues",
    "SJ samples": "Greens", "SJ reads": "Greens",
    "Assembly transcripts": "Oranges",
    "TElocal reads": "Purples",
}

# Below this many shown rows, fill down into the next Screens tier(s) so a
# run whose top Screens value only has a handful of pairs still gets a
# useful table -- see _select_rows(). Not a config key, same reasoning as
# --max-rows.
MIN_ROWS = 50


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


def _select_rows(rows, min_rows, max_rows):
    """Choose which pairs the report renders.

    rows must already be sorted by (-n_screens, gene label, te_id) -- every
    pair sharing a Screens value is therefore a contiguous run. The rule:
    show every pair tied at the top Screens value; if that group alone is
    smaller than min_rows, fill down with whole next-Screens-value tiers
    (never a partial tier) until at least min_rows are shown or the rows
    run out; max_rows is a hard cap that can only ever cut INTO the top
    group itself (a lower tier is only ever added whole or not at all, so
    it never needs a disclosure -- see truncated_top below).

    Returns (selected, top_screens, n_top_group, truncated_top):
      top_screens    the highest n_screens value in the WHOLE file
      n_top_group    how many pairs share top_screens in the WHOLE file,
                     regardless of how many are actually shown
      truncated_top  True only when max_rows had to cut into the top group
                     itself -- the one case where which rows are shown is
                     an alphabetical accident rather than a complete tier,
                     and the section must say so.
    """
    if not rows:
        return [], 0, 0, False

    top_screens = _int(rows[0].get("n_screens"))
    n_top_group = 0
    for r in rows:
        if _int(r.get("n_screens")) != top_screens:
            break
        n_top_group += 1

    if n_top_group > max_rows:
        return rows[:max_rows], top_screens, n_top_group, True

    selected = rows[:n_top_group]
    i = n_top_group
    while len(selected) < min_rows and i < len(rows):
        tier_screens = _int(rows[i].get("n_screens"))
        j = i
        while j < len(rows) and _int(rows[j].get("n_screens")) == tier_screens:
            j += 1
        remaining_capacity = max_rows - len(selected)
        if remaining_capacity <= 0:
            break
        selected += rows[i:j][:remaining_capacity]
        i = j
    return selected, top_screens, n_top_group, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True,
                    help="results/chimera/candidates.tsv.gz")
    ap.add_argument("--gene-names", default=None,
                    help="gene_id_to_name.tsv.gz, to label rows by symbol")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-rows", type=int, default=500,
                    help="Hard safety cap on rendered rows -- not the "
                    "normal row count, see _select_rows()'s docstring.")
    ap.add_argument("--source-path", default="results/chimera/candidates.tsv.gz")
    ap.add_argument("--explorer-path",
                    default="results/chimera/candidates_explorer.html",
                    help="chimera_candidates_explorer.R's output -- the full, "
                    "sortable/filterable table, no top-N cap")
    args = ap.parse_args()

    rows = list(load(args.evidence))
    symbols = load_symbols(args.gene_names)
    for r in rows:
        r["_gene_label"] = symbols.get(r.get("gene_id", "."), r.get("gene_id", "."))

    # Selection: most screens first, ties broken alphabetically -- a sort
    # that asserts nothing beyond screen count. candidates.tsv.gz's own
    # order additionally breaks ties by n_screen_evidence/n_corroboration
    # (see chimera_evidence.py), but this table no longer shows either of
    # those as columns, so re-using the file's finer tie-break here would
    # silently choose which rows are shown on signals the reader can no
    # longer see or verify -- the report and the file therefore do NOT
    # always agree on row order within a tied Screens value any more.
    # candidates.tsv.gz remains the source of truth for the full, finely
    # sorted catalogue.
    rows.sort(key=lambda r: (
        -_int(r.get("n_screens")), r["_gene_label"], r.get("te_id", "."),
    ))

    top, top_screens, n_top_group, truncated_top = _select_rows(
        rows, MIN_ROWS, args.max_rows
    )

    data = {}
    for r in top:
        gene = r["_gene_label"]
        # " | ", never " / ": MultiQC cleans table row names like filenames
        # and splits on "/", taking the basename -- which silently dropped the
        # gene and left only the TE id. Measured; guard 50 pins it.
        key = f"{gene} | {r.get('te_id', '.')}"
        # a duplicate key would silently drop a row; disambiguate with gene_id
        if key in data:
            key = f"{key} ({r.get('gene_id', '.')})"
        entry = {}
        for col, header, _desc, kind in COLUMNS:
            if col == "gene_label":
                value = gene
            elif col == "replicated":
                # Mirrors chimera_evidence.py's own multi_sample rule
                # exactly (max(cr_max_samples, sj_max_samples) > 1) rather
                # than recomputing it from the raw sample counts here.
                value = ("yes" if "multi_sample"
                         in r.get("corroboration", ".").split(",") else "no")
            elif col == "types":
                calls = {t for c in ("cr_chimera_types", "sj_chimera_types",
                                     "assembly_chimera_types")
                         for t in r.get(c, ".").split(",")
                         if t in CHIMERA_CALL_TYPES}
                value = ",".join(sorted(calls)) or "."
            else:
                value = r.get(col, ".")
            if kind == "intna":
                # Leave the cell UNSET when the measurement was not taken.
                # _int(".") is 0, which would read as "TElocal ran and found
                # no reads" -- a different statement from "TElocal did not
                # run". MultiQC renders a missing key as an empty cell.
                if str(value) not in (".", ""):
                    entry[header] = _int(value)
            else:
                entry[header] = _int(value) if kind == "int" else str(value)
        data[key] = entry

    # placement guarantees column order regardless of dict/JSON key-order
    # quirks -- explicit rather than relying on insertion order alone.
    headers = {}
    for i, (_col, header, desc, kind) in enumerate(COLUMNS):
        h = {"title": header, "description": desc, "placement": i * 10}
        if kind in ("int", "intna"):
            h["format"] = "{:,.0f}"
            h["min"] = 0
        if header in HIDDEN_HEADERS:
            h["hidden"] = True
        if header in SCALE_BY_HEADER:
            h["scale"] = SCALE_BY_HEADER[header]
        headers[header] = h

    # Only the top-group-truncation case involves any alphabetical choice
    # (see _select_rows()'s docstring) -- every other case shows either the
    # whole top group, or the whole top group plus whole lower tiers, so it
    # is described plainly instead.
    if truncated_top:
        tie_note = (
            f"{n_top_group:,} pairs share the top Screens value "
            f"({top_screens}) -- only the alphabetically first "
            f"{args.max_rows:,} of those are shown here; the rest are in "
            "the explorer. "
        )
    else:
        extra = len(top) - n_top_group
        if extra > 0:
            tie_note = (
                f"All <strong>{n_top_group:,}</strong> pairs found by the "
                f"most screens (Screens = {top_screens}) are shown, plus "
                f"the next <strong>{extra:,}</strong> pairs by Screens "
                "value, so the table stays a useful size. "
            )
        elif rows:
            tie_note = (
                f"All <strong>{n_top_group:,}</strong> pairs found by the "
                f"most screens (Screens = {top_screens}) are shown. "
            )
        else:
            tie_note = ""

    # Kept deliberately short: the full argument lives in "How to weigh this
    # evidence" one section down (see report_section_order in
    # multiqc_config.yaml: Candidates renders first), and repeating it here
    # buried the table. What
    # must survive any trim is the disclosure that the default order is an
    # unweighted count rather than a ranking -- guard 50 pins the phrases
    # "not a score", "sort" and "validate candidates manually".
    note = (
        "<p>The <strong>{n_shown:,}</strong> of {n_total:,} gene-TE pairs "
        "shown here, from <code>{src}</code>: sorted by <strong>Screens"
        "</strong> descending, ties broken alphabetically by gene then TE. "
        "{tie_note}"
        "<strong>Click any column header to sort</strong> -- every column "
        "is an individual, unweighted signal, <strong>not a score</strong>; "
        "see <strong>How to weigh this evidence</strong> below, and "
        "validate candidates manually. For all <strong>{n_total:,}</strong> "
        "pairs, filterable, no MultiQC needed: <code>{explorer}</code>.</p>"
    ).format(n_shown=len(top), n_total=len(rows), src=args.source_path,
             tie_note=tie_note, explorer=args.explorer_path)
    note_help = (
        "<p>Ordered by <strong>Screens</strong> only -- the number of "
        "independent screens (cr / assembly / sj) that found a pair -- "
        "then alphabetically by gene, then TE. Screens is a count, not a "
        "score: every other signal (CR motif, SJ motif, Assembly strand, "
        "Replicated, TElocal reads) is its own column instead of being "
        "folded into a second count, so nothing beyond screen agreement "
        "decides which rows are shown first.</p>"
        "<p>Which rows: every pair tied at the top Screens value is shown, "
        "in full. If that group is smaller than a useful table, whole next-"
        "Screens-value tiers are added until it is. A hard cap "
        "(--max-rows) only ever bites when the top group itself is bigger "
        "than the cap -- the one case where which rows are shown is an "
        "alphabetical accident, and the section says so explicitly when it "
        "happens.</p>"
    )

    if top:
        body = {
            "plot_type": "table",
            "pconfig": {
                "id": "chimera_candidates_table_plot",
                "title": "Gene-TE chimera candidates",
                "col1_header": "Gene | TE",
                # defaultsort is what actually sets the opening order.
                # sort_rows: False does NOT survive -- MultiQC re-populates its
                # camelCase alias sortRows from the default (True), so the rows
                # arrive alphabetised by name whatever this says. Measured.
                # Stating the intended sort explicitly is the reliable route.
                "sort_rows": False,
                "defaultsort": [
                    {"column": "Screens", "direction": "desc"},
                ],
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
        "helptext": note_help,
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
