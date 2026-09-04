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
rungs, and the pipeline's own measurements contradicted its top rung (see
chimera_evidence_heatmap.py: cross-screen agreement sits near its chance
rate).

Two columns summarise what was observed, both deliberately unweighted:

  evidence    the set of evidence types present, comma-joined ("." if none).
              Unordered and unweighted -- the flags are names, not points:

                canonical              a recognised splice motif on at least
                                       one junction
                multi_sample           seen in more than one sample
                both_screens           called by BOTH screens
                assembly_strand_match  the assembled transcript's strand
                                       agrees with the gene's
                telocal_expressed      TElocal reports the TE locus as
                                       expressed (unresolved signal -- see
                                       below)

  n_evidence  how many of those flags are set.  A COUNT OF EVIDENCE TYPES,
              NOT A CONFIDENCE SCORE.  It weights every flag equally for the
              specific reason that no weighting has been validated here.  The
              file is sorted by it only so the order is deterministic and the
              densely-evidenced rows are easy to find; it is not a claim that
              those rows are correct.

              Known bias, stated because the number invites over-reading:
              n_evidence structurally favours pairs the assembly screen found,
              since both_screens and assembly_strand_match are unreachable
              without assembly support.  It is a tally of what was observed,
              not a comparison of candidates.

Read depth (junction_reads / junction_events) is deliberately NOT an
evidence flag, because it looks like support and is not: the metric most
inflated by artifacts -- a hot PCR chimera is often the deepest event in a
run. It is still reported as a column.

TElocal expression of the TE locus (telocal_expressed) IS counted, but its
standing is unresolved, not confirmed. Measured on a real 4-sample mouse
run it looked like the opposite of support: 91% of junction-side pairs had
an expressed locus, so it discriminated nothing, and the canonical rate was
LOWER where the TE was expressed (6.7% at telocal_count > 10 vs 10.2% at
<= 10, n = 19,503 events) -- mechanistically unsurprising (a highly
expressed locus yields more reads and so more chances for template
switching), but not evidence of a real chimera either. One small run isn't
enough to demote a signal on, so it stays a flag until that correlation is
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
    "found_by", "evidence", "n_evidence",
    "junction_events", "junction_reads", "junction_max_samples",
    "junction_canonical", "junction_chimera_types",
    "telocal_active", "telocal_count",
    "assembly_transcripts", "assembly_chimera_types",
    "assembly_strand_match", "assembly_transcript_ids",
]


def _blank():
    return {
        "te_subfamily": ".", "te_family": ".", "te_class": ".",
        "junction_events": 0, "junction_reads": 0, "junction_max_samples": 0,
        "junction_canonical": "no", "junction_types": set(),
        "telocal_active": ".", "telocal_count": 0,
        "assembly_transcripts": 0, "assembly_types": set(),
        "assembly_strand_match": ".", "assembly_tids": [],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junction", required=True,
                    help="results/chimera/reads/te-gene-chimeras.tsv.gz")
    ap.add_argument("--assembly", default=None,
                    help="results/chimera/assembly/transcripts.tsv.gz "
                         "(omit when the assembly screen is disabled)")
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
        p["junction_events"] += 1
        p["junction_reads"] += _int(r.get("total_reads"))
        p["junction_max_samples"] = max(
            p["junction_max_samples"], _int(r.get("n_samples"))
        )
        if r.get("canonical") == "yes":
            p["junction_canonical"] = "yes"
        if r.get("chimera_type", ".") != ".":
            p["junction_types"].add(r["chimera_type"])
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
        # derived from. MAX, not sum: telocal_count is one locus's count in one
        # sample, and the same locus recurs across a pair's events, so summing
        # would multiply one measurement by how many junctions happened to hit
        # it. Max reads as "the most this locus was expressed in any sample".
        p["telocal_count"] = max(p["telocal_count"], _int(r.get("telocal_count", 0)))

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

    rows = []
    for (gene, te), p in pairs.items():
        in_junction = p["junction_events"] > 0
        in_assembly = p["assembly_transcripts"] > 0
        found_by = (
            "both" if in_junction and in_assembly
            else "reads" if in_junction
            else "assembly"
        )
        # Names, not points. Order here is presentational only -- nothing
        # downstream may treat position in this list as a weight.
        flags = []
        if p["junction_canonical"] == "yes":
            flags.append("canonical")
        if p["junction_max_samples"] > 1:
            flags.append("multi_sample")
        if found_by == "both":
            flags.append("both_screens")
        if p["assembly_strand_match"] == "yes":
            flags.append("assembly_strand_match")
        # TElocal expression COUNTS as evidence, deliberately. One 4-sample
        # run suggested it discriminates nothing (see the evidence guide), but
        # a single small experiment is not enough to demote a signal: the
        # correlation between junction-side pairs and locus expression has not
        # been tested properly yet. It stays a flag until it has been.
        if p["telocal_active"] == "yes":
            flags.append("telocal_expressed")

        rows.append({
            "gene_id": gene, "te_id": te,
            "te_subfamily": p["te_subfamily"], "te_family": p["te_family"],
            "te_class": p["te_class"],
            "found_by": found_by,
            "evidence": ",".join(flags) or ".",
            "n_evidence": len(flags),
            "junction_events": p["junction_events"],
            "junction_reads": p["junction_reads"],
            "junction_max_samples": p["junction_max_samples"],
            "junction_canonical": p["junction_canonical"],
            "junction_chimera_types": ",".join(sorted(p["junction_types"])) or ".",
            "telocal_active": p["telocal_active"],
            # "." rather than 0 when TElocal never ran, so "not measured" stays
            # distinguishable from "measured, no reads" -- the same distinction
            # telocal_active already draws.
            "telocal_count": ("." if p["telocal_active"] == "."
                              else p["telocal_count"]),
            "assembly_transcripts": p["assembly_transcripts"],
            "assembly_chimera_types": ",".join(sorted(p["assembly_types"])) or ".",
            "assembly_strand_match": p["assembly_strand_match"],
            "assembly_transcript_ids": ",".join(p["assembly_tids"]) or ".",
        })

    # Deterministic order: densest evidence first, then alphabetical. This is
    # a sort, not a verdict -- n_evidence counts flags without weighting them,
    # and gene/te break ties so two runs of the same data produce byte-identical
    # files. Nothing here says a high-n_evidence pair is real.
    rows.sort(key=lambda r: (-r["n_evidence"], r["gene_id"], r["te_id"]))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(OUT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in OUT_COLUMNS) + "\n")

    composition = {}
    for r in rows:
        for flag in r["evidence"].split(","):
            if flag != ".":
                composition[flag] = composition.get(flag, 0) + 1
    summary = ", ".join(
        f"{flag}: {n}" for flag, n in sorted(composition.items())
    )
    print(f"chimera evidence: {len(rows)} gene-TE pairs ({summary or 'no evidence flags'}) "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
