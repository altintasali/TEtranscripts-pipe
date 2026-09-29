#!/usr/bin/env python3
"""One unified gene-TE chimera table: every line of evidence for a candidate
in one row.

The problem this solves: the two screens key their output differently -- the
junction screen is BREAKPOINT-keyed (one row per chimeric junction, many per
locus) and the assembly screen is TRANSCRIPT-keyed (one row per assembled
transcript) -- so neither table can be read as "the candidate list", and a
user comparing them by eye is doing a join by hand.

The one key both screens genuinely share is the **(gene_id, te_id) pair**,
where te_id is the individual TE insertion (transcript_id, e.g. L1PA2_dup1)
rather than the subfamily.  That is what this collapses to: one row per
gene-TE pair, with the per-screen detail aggregated into it.  It is a
strictly coarser view than either source table -- to see the individual
breakpoints or transcripts behind a row, go back to
te-gene-chimeras.tsv.gz / candidates.tsv.gz using the same pair.

This table reports evidence.  It does NOT score or rank candidates, and it
does not decide which are real -- that is a manual call, made against these
columns.  An earlier version carried a four-tier confidence ladder; it was
removed because no experiment here established the relative weight of its
rungs, and an early project measurement contradicted its top rung (see
chimera_evidence_guide_mqc.py's "Screens / Found by" row, and
docs/chimera-evidence.md for the measurement itself, its cohort and its
restrictions).

Six columns summarise what was observed, all deliberately unweighted, in
three pairs -- which screens found it / how many, what SCREEN-BOUND quality
signal each of those screens gave it / how many, and what CROSS-CUTTING
corroboration exists / how many:

  found_by          which of the three independent screens (cr / assembly /
                     sj) contributed a row for this pair, "+"-joined (e.g.
                     "cr+sj"). Not itself weighted or ordered.

  n_screens         how many of those three screens found the pair (1-3). A
                     plain count derived from found_by -- NOT summed into
                     either count below. It used to be (see git history):
                     two flags, both_screens and all_three_screens,
                     duplicated this same found-by-N-screens fact and were
                     counted a second time, so a pair found by all 3 screens
                     got +2 just from that overlap on top of everything
                     else. n_screens replaces both, without double-counting.

  screen_evidence    quality signals tied to a SPECIFIC screen, comma-joined
                     ("." if none) -- each one can only be set if that
                     screen's own rows are present for this pair, so this
                     set is always a subset of what found_by names:

                       cr_canonical           a recognised splice motif on
                                              at least one chimeric-junction
                                              read (cr screen)
                       sj_canonical           a recognised splice motif on
                                              at least one SJ.out.tab
                                              junction (sj screen) -- kept
                                              separate from cr_canonical:
                                              structurally independent
                                              measurements
                       assembly_strand_match  the assembled transcript's
                                              strand agrees with the gene's
                                              (assembly screen)

  n_screen_evidence  how many of those three are set (0-3). Because each one
                     requires its own screen to have found the pair,
                     n_screen_evidence <= n_screens ALWAYS holds -- unlike
                     corroboration below, this count really is bounded by
                     how many screens fired.

                     Residual bias, stated because the number invites
                     over-reading: n_screen_evidence still favours pairs the
                     assembly screen found, since assembly_strand_match is
                     unreachable without assembly support. It is a tally of
                     what was observed, not a comparison of candidates.

  corroboration      signals that do NOT belong to any one screen,
                     comma-joined ("." if none) -- unlike screen_evidence,
                     these can be present even for a pair found by only one
                     screen, which is exactly why they are kept in a
                     separate column rather than folded into the same
                     count as the screen-bound flags above:

                       multi_sample        seen in more than one sample --
                                           can fire from a SINGLE screen
                                           alone, no second screen required
                       telocal_expressed   TElocal reports the TE locus as
                                           expressed (unresolved signal --
                                           see below) -- from TElocal, a
                                           FOURTH data source that is not
                                           one of the three detection
                                           screens at all

  n_corroboration    how many of those two are set (0-2). Deliberately NOT
                     comparable to n_screens or n_screen_evidence -- a pair
                     found by exactly one screen can still reach
                     n_corroboration == 2.

Neither count is a confidence score. Both weight every flag inside them
equally, for the specific reason that no weighting has been validated here.
The file is sorted by (n_screens, n_screen_evidence, n_corroboration), all
descending, only so the order is deterministic and the densely-evidenced
rows are easy to find; it is not a claim that those rows are correct.

Read depth (cr_reads / cr_events) is deliberately NOT a
flag in either count, because it looks like support and is not: the metric
most inflated by artifacts -- a hot PCR chimera is often the deepest event
in a run. It is still reported as a column.

cr_gene_te_distance is reported, never counted: the genomic gap between
the gene and the TE for pairs the chimeric-reads screen found -- "trans"
(different chromosomes), 0 (the TE overlaps the gene's span) or bp; "." for
pairs that screen did not find. On a real run most chimeric-read gene<->TE
pairs were trans or >200 kb apart, the random-partner pattern of template
switching / chimeric ligation rather than a TE driving that gene; this
column lets them be told apart from local events.

The three *_chimera_types columns can include antisense_to_gene: a call
whose transcript runs on the strand opposite its gene (see
chimera_exon_context.ANTISENSE_TO_GENE).

TElocal expression of the TE locus (telocal_expressed) IS counted, but its
standing is not validated, not confirmed. An early project measurement
(see docs/chimera-evidence.md) looked like the opposite of support rather
than for it -- mechanistically unsurprising either way (a highly expressed
locus yields more reads and so more chances for template switching, but
also more chances to actually observe a real chimera), and not enough on
its own to demote a signal, so it stays a flag until that correlation is
tested properly across more data; see chimera_evidence_guide_mqc.py for the
report-facing version of this caveat. telocal_count/telocal_active are
always reported regardless of the flag.

Relative weight of the flags is exactly what has NOT been established, so
the report states what is known about each one (see
chimera_evidence_guide_mqc.py) rather than combining them.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def load(path):
    """Rows of a TSV as dicts, keyed by header name."""
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


OUT_COLUMNS = [
    "gene_id", "te_id", "te_subfamily", "te_family", "te_class",
    "found_by", "n_screens",
    "screen_evidence", "n_screen_evidence",
    "corroboration", "n_corroboration",
    "cr_events", "cr_reads", "cr_max_samples",
    "cr_canonical", "cr_chimera_types", "cr_gene_te_distance",
    "telocal_active", "telocal_count", "telocal_locus",
    "assembly_transcripts", "assembly_chimera_types",
    "assembly_strand_match", "assembly_transcript_ids",
    "sj_events", "sj_reads", "sj_max_samples",
    "sj_canonical", "sj_chimera_types",
]


def _blank():
    return {
        "te_subfamily": ".", "te_family": ".", "te_class": ".",
        "cr_events": 0, "cr_reads": 0, "cr_max_samples": 0,
        "cr_canonical": "no", "junction_types": set(),
        "cr_gene_te_distance": ".",
        "telocal_active": ".", "telocal_count": 0, "telocal_locus": ".",
        "assembly_transcripts": 0, "assembly_types": set(),
        "assembly_strand_match": ".", "assembly_tids": [],
        "sj_events": 0, "sj_reads": 0, "sj_max_samples": 0,
        "sj_canonical": "no", "sj_types": set(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junction", required=True,
                    help="results/chimera/chimeric_reads/te-gene-chimeras.tsv.gz")
    ap.add_argument("--assembly", default=None,
                    help="results/chimera/assembly/transcripts.tsv.gz "
                         "(omit when the assembly screen is disabled)")
    ap.add_argument("--sj", default=None,
                    help="results/chimera/splice_junctions/te-gene-junctions.tsv.gz "
                         "(omit when the SJ.out.tab screen is disabled)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pairs = {}

    for r in load(args.junction):
        gene, te = r.get("gene_id", "."), r.get("te_id", ".")
        if gene in (".", "") or te in (".", ""):
            continue
        p = pairs.setdefault((gene, te), _blank())
        # Annotation is per-insertion, so every row for a pair agrees; take
        # the first non-"." rather than overwriting with a later "."
        for col in ("te_subfamily", "te_family", "te_class"):
            if p[col] == "." and r.get(col, ".") != ".":
                p[col] = r[col]
        p["cr_events"] += 1
        p["cr_reads"] += _int(r.get("total_reads"))
        p["cr_max_samples"] = max(
            p["cr_max_samples"], _int(r.get("n_samples"))
        )
        if r.get("canonical") == "yes":
            p["cr_canonical"] = "yes"
        if r.get("chimera_type", ".") != ".":
            p["junction_types"].add(r["chimera_type"])
        # gene<->TE distance is a property of the pair (same gene, same TE
        # insertion), so every event agrees; "." for tables written before
        # the column existed.
        if p["cr_gene_te_distance"] == ".":
            p["cr_gene_te_distance"] = r.get("gene_te_distance", ".") or "."
        # "yes" for the pair if ANY event's TE locus is called expressed;
        # stays "." (not "no") when telocal never ran, so an absent check is
        # distinguishable from a negative one.
        active = r.get("telocal_active", ".")
        if active == "yes":
            p["telocal_active"] = "yes"
        elif active == "no" and p["telocal_active"] == ".":
            p["telocal_active"] = "no"
        # chimera_telocal_annotate.py already measured the locus's read count
        # and this step used to throw it away, keeping only the boolean it was
        # derived from. MAX, not sum: telocal_count is already a COHORT-WIDE
        # total for that locus (TelocalIndex.build in chimera_telocal_index.py
        # sums every sample's cntTable into one shared index before any
        # per-sample annotate job ever sees it -- it is NOT a per-sample
        # figure despite the name), and the same locus recurs across a pair's
        # events, so summing here would multiply that one cohort total by how
        # many junctions happened to hit it. Max just re-selects the same
        # constant value rather than inflating it.
        p["telocal_count"] = max(p["telocal_count"], _int(r.get("telocal_count", 0)))
        # First non-"." locus key seen, same convention as te_subfamily/
        # te_family/te_class above -- it's the TElocal cntTable key for this
        # TE copy, breakpoint-deterministic like the rest of the annotation,
        # so every row for a pair agrees. Needed to join a real cohort-total
        # TElocal count from results/telocal/counts_matrix.tsv.gz downstream
        # (candidates_explorer.html) -- which will agree with telocal_count
        # above, since both are cohort-wide sums computed two different ways;
        # the matrix join is kept as the one actually displayed there.
        if p["telocal_locus"] == "." and r.get("telocal_locus", ".") != ".":
            p["telocal_locus"] = r["telocal_locus"]

    if args.assembly:
        for r in load(args.assembly):
            gene, te = r.get("matched_gene_id", "."), r.get("te_id", ".")
            if gene in (".", "") or te in (".", ""):
                continue
            p = pairs.setdefault((gene, te), _blank())
            for col in ("te_subfamily", "te_family", "te_class"):
                if p[col] == "." and r.get(col, ".") != ".":
                    p[col] = r[col]
            p["assembly_transcripts"] += 1
            if r.get("chimera_type", ".") != ".":
                p["assembly_types"].add(r["chimera_type"])
            if r.get("strand_match") == "yes":
                p["assembly_strand_match"] = "yes"
            elif p["assembly_strand_match"] == ".":
                p["assembly_strand_match"] = r.get("strand_match", ".")
            p["assembly_tids"].append(r.get("transcript_id", "."))

    if args.sj:
        for r in load(args.sj):
            gene, te = r.get("gene_id", "."), r.get("te_id", ".")
            if gene in (".", "") or te in (".", ""):
                continue
            p = pairs.setdefault((gene, te), _blank())
            for col in ("te_subfamily", "te_family", "te_class"):
                if p[col] == "." and r.get(col, ".") != ".":
                    p[col] = r[col]
            p["sj_events"] += 1
            p["sj_reads"] += _int(r.get("total_reads"))
            p["sj_max_samples"] = max(
                p["sj_max_samples"], _int(r.get("n_samples"))
            )
            if r.get("canonical") == "yes":
                p["sj_canonical"] = "yes"
            if r.get("chimera_type", ".") != ".":
                p["sj_types"].add(r["chimera_type"])

    rows = []
    for (gene, te), p in pairs.items():
        in_junction = p["cr_events"] > 0
        in_assembly = p["assembly_transcripts"] > 0
        in_sj = p["sj_events"] > 0
        sources = []
        if in_junction:
            sources.append("cr")
        if in_assembly:
            sources.append("assembly")
        if in_sj:
            sources.append("sj")
        # "+"-joined screen tags, matching the same abbreviations used in
        # the evidence flags (cr_canonical / sj_canonical). No special-cased
        # word for any one combination: an earlier version wrote "both" for
        # exactly {reads, assembly} (back when those were the only two
        # screens), but keeping one particular combination as a word while
        # every other combination is "+"-joined was its own inconsistency
        # once a third screen (sj) existed -- and "both" doesn't even parse
        # once the reads screen's tag is "cr" instead of "reads".
        found_by = "+".join(sources)
        # How many of the three screens found this pair -- NOT summed into
        # n_evidence below (see the module docstring: this replaces the old
        # both_screens/all_three_screens flags, which duplicated this same
        # fact and got double-counted).
        n_screens = int(in_junction) + int(in_assembly) + int(in_sj)
        # Names, not points. Order here is presentational only -- nothing
        # downstream may treat position in this list as a weight.
        #
        # Split into two independent counts because they answer different
        # questions and are NOT nested: screen_evidence flags can only fire
        # if their own screen found the pair (so n_screen_evidence <=
        # n_screens always), but corroboration flags can fire regardless of
        # how many screens found it -- folding them into one count made
        # n_evidence look like it should relate to n_screens when it
        # structurally could not (a single-screen pair could still out-count
        # a three-screen one).
        screen_evidence = []
        if p["cr_canonical"] == "yes":
            screen_evidence.append("cr_canonical")
        # Kept separate from "cr_canonical" (reads-screen junction type)
        # rather than merged: they are two structurally independent
        # measurements (STAR chimeric-junction typing vs. STAR SJ.out.tab
        # motif), and folding them into one flag would hide which one
        # actually fired -- against this file's own reason for having a
        # "source" column at all.
        if p["sj_canonical"] == "yes":
            screen_evidence.append("sj_canonical")
        if p["assembly_strand_match"] == "yes":
            screen_evidence.append("assembly_strand_match")

        corroboration = []
        if max(p["cr_max_samples"], p["sj_max_samples"]) > 1:
            corroboration.append("multi_sample")
        # TElocal expression COUNTS as corroboration, deliberately. One
        # 4-sample run suggested it discriminates nothing (see the evidence
        # guide), but a single small experiment is not enough to demote a
        # signal: the correlation between junction-side pairs and locus
        # expression has not been tested properly yet. It stays a flag
        # until it has been.
        if p["telocal_active"] == "yes":
            corroboration.append("telocal_expressed")

        rows.append({
            "gene_id": gene, "te_id": te,
            "te_subfamily": p["te_subfamily"], "te_family": p["te_family"],
            "te_class": p["te_class"],
            "found_by": found_by,
            "n_screens": n_screens,
            "screen_evidence": ",".join(screen_evidence) or ".",
            "n_screen_evidence": len(screen_evidence),
            "corroboration": ",".join(corroboration) or ".",
            "n_corroboration": len(corroboration),
            "cr_events": p["cr_events"],
            "cr_reads": p["cr_reads"],
            "cr_max_samples": p["cr_max_samples"],
            "cr_canonical": p["cr_canonical"],
            "cr_chimera_types": ",".join(sorted(p["junction_types"])) or ".",
            "cr_gene_te_distance": p["cr_gene_te_distance"],
            "telocal_active": p["telocal_active"],
            # "." rather than 0 when TElocal never ran, so "not measured" stays
            # distinguishable from "measured, no reads" -- the same distinction
            # telocal_active already draws.
            "telocal_count": ("." if p["telocal_active"] == "."
                              else p["telocal_count"]),
            "telocal_locus": p["telocal_locus"],
            "assembly_transcripts": p["assembly_transcripts"],
            "assembly_chimera_types": ",".join(sorted(p["assembly_types"])) or ".",
            "assembly_strand_match": p["assembly_strand_match"],
            "assembly_transcript_ids": ",".join(p["assembly_tids"]) or ".",
            "sj_events": p["sj_events"],
            "sj_reads": p["sj_reads"],
            "sj_max_samples": p["sj_max_samples"],
            "sj_canonical": p["sj_canonical"],
            "sj_chimera_types": ",".join(sorted(p["sj_types"])) or ".",
        })

    # Deterministic order: most screens first, then densest screen-bound
    # evidence, then most corroboration, then alphabetical. This is a sort,
    # not a verdict -- each count is flags without weighting them, and
    # gene/te break ties so two runs of the same data produce byte-identical
    # files. Nothing here says a row sorted first is real.
    rows.sort(key=lambda r: (
        -r["n_screens"], -r["n_screen_evidence"], -r["n_corroboration"],
        r["gene_id"], r["te_id"],
    ))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(OUT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in OUT_COLUMNS) + "\n")

    composition = {}
    for r in rows:
        for flag in r["screen_evidence"].split(",") + r["corroboration"].split(","):
            if flag != ".":
                composition[flag] = composition.get(flag, 0) + 1
    summary = ", ".join(
        f"{flag}: {n}" for flag, n in sorted(composition.items())
    )
    print(f"chimera evidence: {len(rows)} gene-TE pairs ({summary or 'no evidence flags'}) "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
