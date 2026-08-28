#!/usr/bin/env python3
"""One unified gene-TE chimera table: every line of evidence for a candidate
in one row, with a confidence tier.

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

Confidence tiers mirror the ladder documented in the wiki's "Interpreting
Your Results" page, and are the MAXIMUM level a pair reaches:

  1  a gene-TE call exists at all (the floor for any row here)
  2  + seen in more than one sample
  3  + a recognised splice motif on at least one junction (canonical)
  4  + called by BOTH screens (methods with opposite blind spots agree)

Two things are deliberately NOT on the ladder, both because they look like
support and are not:

Read depth.  It is the metric most inflated by artifacts -- a hot PCR
chimera is often the deepest event in a run.

TElocal expression of the TE locus.  This WAS a tier, on the assumption
that independent evidence the TE is transcribed corroborates the chimera.
Measured on a real 4-sample mouse run it does the opposite: 91% of
junction-side pairs had an expressed locus, so it discriminated nothing,
and the canonical rate was LOWER where the TE was expressed (6.7% at
telocal_count > 10 vs 10.2% at <= 10, n = 19,503 events).  That is
mechanistically unsurprising -- a highly expressed locus yields more reads
and so more chances for template switching -- but it means expression is
context, not support.  telocal_count/telocal_active are still reported;
they just no longer promote anything.

Ordering note: canonical outranks reproducibility because the dominant
artifact mode here is sequence-driven template switching, which recurs
across libraries -- so being seen in several samples does not rule it out,
while a splice motif does.
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
    "found_by", "confidence_tier",
    "junction_events", "junction_reads", "junction_max_samples",
    "junction_canonical", "junction_chimera_types", "telocal_active",
    "assembly_transcripts", "assembly_chimera_types",
    "assembly_strand_match", "assembly_transcript_ids",
]


def _blank():
    return {
        "te_subfamily": ".", "te_family": ".", "te_class": ".",
        "junction_events": 0, "junction_reads": 0, "junction_max_samples": 0,
        "junction_canonical": "no", "junction_types": set(),
        "telocal_active": ".",
        "assembly_transcripts": 0, "assembly_types": set(),
        "assembly_strand_match": ".", "assembly_tids": [],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junction", required=True,
                    help="results/chimera_junction/te-gene-chimeras.tsv.gz")
    ap.add_argument("--assembly", default=None,
                    help="results/chimera_assembly/candidates.tsv.gz "
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
            else "junction" if in_junction
            else "assembly"
        )
        tier = 1
        if p["junction_max_samples"] > 1:
            tier = max(tier, 2)
        if p["junction_canonical"] == "yes":
            tier = max(tier, 3)
        if found_by == "both":
            tier = max(tier, 4)

        rows.append({
            "gene_id": gene, "te_id": te,
            "te_subfamily": p["te_subfamily"], "te_family": p["te_family"],
            "te_class": p["te_class"],
            "found_by": found_by, "confidence_tier": tier,
            "junction_events": p["junction_events"],
            "junction_reads": p["junction_reads"],
            "junction_max_samples": p["junction_max_samples"],
            "junction_canonical": p["junction_canonical"],
            "junction_chimera_types": ",".join(sorted(p["junction_types"])) or ".",
            "telocal_active": p["telocal_active"],
            "assembly_transcripts": p["assembly_transcripts"],
            "assembly_chimera_types": ",".join(sorted(p["assembly_types"])) or ".",
            "assembly_strand_match": p["assembly_strand_match"],
            "assembly_transcript_ids": ",".join(p["assembly_tids"]) or ".",
        })

    # Best candidates first, so head-ing the file is a useful triage step.
    # Reads break ties last, matching the ladder's deliberate ordering.
    rows.sort(
        key=lambda r: (
            r["confidence_tier"], r["junction_max_samples"],
            r["assembly_transcripts"], r["junction_reads"],
        ),
        reverse=True,
    )

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(OUT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in OUT_COLUMNS) + "\n")

    by_tier = {}
    for r in rows:
        by_tier[r["confidence_tier"]] = by_tier.get(r["confidence_tier"], 0) + 1
    summary = ", ".join(f"tier {t}: {n}" for t, n in sorted(by_tier.items(), reverse=True))
    print(f"chimera evidence: {len(rows)} gene-TE pairs ({summary or 'none'}) "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
