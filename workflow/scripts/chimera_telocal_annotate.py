#!/usr/bin/env python3
"""Cross-reference chimera junction tables with TElocal locus-level expression.

Strategy A: For each gene<->TE chimera event, check whether the TE breakpoint
overlaps a TElocal locus that has nonzero read counts.  Add three columns:
  telocal_count    read count from the best-matching TElocal locus (0 if none)
  telocal_locus    the TElocal locus key, or "." if no overlap
  telocal_active   "yes" if telocal_count > 0, else "no"

Strategy B: When multiple TEs overlap a chimera breakpoint (visible in the
donor_hits / acceptor_hits columns), use TElocal expression to pick the
transcribed TE instead of the alphabetically-first one.  If the swap happens,
te_id is updated and te_refined_by_telocal is set to "yes".

Non-chimera events (direction not gene_to_te / te_to_gene) receive the default
values (0 / "." / "no" / ".").

Usage:
  chimera_telocal_annotate.py \\
      --junctions {sample}_junctions.tsv.gz \\
      --telocal-index telocal_index.pkl.gz \\
      --sample-name {sample} \\
      --breakpoint-tolerance 0 \\
      --out {sample}_junctions_with-telocal.tsv.gz \\
      --te-out {sample}_junctions_te-gene-chimeras_with-telocal.tsv.gz

--out keeps every junction (same rows as --junctions, plus the four columns
above).  --te-out is the optional gene<->TE subset of that same table --
same filter classify_chimera_junctions.py's own --te-out applies, so the four
per-sample tables are the two row-sets (all / gene-TE) x the two column-sets
(without / with TElocal).

The interval index (--telocal-index) is built ONCE for the whole cohort by
build_chimera_telocal_index.py (chimera_telocal_index rule) -- see that
script and chimera_telocal_index.py for why: every per-sample job used to
re-parse and re-sum all samples' cntTables from scratch, an Nx redundant
rebuild of the identical index.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chimera_telocal_index import TelocalIndex, telocal_te_id
from gz_io import open_read, open_write


def parse_hits_csv(csv_str):
    """Parse 'gene:X,Y|te:A,B' hit string to list of TE IDs."""
    if not csv_str or csv_str == ".":
        return []
    te_part = csv_str.split("|te:")[-1] if "|te:" in csv_str else csv_str
    if te_part == "" or te_part == ".":
        return []
    return [x for x in te_part.split(",") if x]


def best_telocal_hit(hits, te_family=None):
    """Pick the best TElocal locus from overlapping hits.

    Preference: matching te_family, then highest count.
    Returns (locus_key, count) or (None, 0) if empty.
    """
    if not hits:
        return None, 0
    if te_family:
        fam_matches = [h for h in hits if h[3] == te_family]
        if fam_matches:
            hits = fam_matches
    best = max(hits, key=lambda h: h[5])
    return best[2], best[5]


def annotate_event(row, telocal_index, breakpoint_tolerance):
    """Add telocal columns to a single junction row (dict). Returns the dict."""
    direction = row.get("direction", "")
    donor_chrom = row.get("donor_chrom", ".")
    acceptor_chrom = row.get("acceptor_chrom", ".")

    row["telocal_count"] = "0"
    row["telocal_locus"] = "."
    row["telocal_active"] = "no"
    row["te_refined_by_telocal"] = "no"

    if direction not in ("gene_to_te", "te_to_gene"):
        return row

    if direction == "gene_to_te":
        te_chrom = acceptor_chrom
        te_bp = int(row.get("acceptor_breakpoint", 0))
    else:
        te_chrom = donor_chrom
        te_bp = int(row.get("donor_breakpoint", 0))

    tol = breakpoint_tolerance
    q_start = te_bp - 1 - tol
    q_end = te_bp + 1 + tol

    hits = telocal_index.overlapping(te_chrom, q_start, q_end)
    te_family = row.get("te_family", ".")
    locus_key, count = best_telocal_hit(hits, te_family)

    row["telocal_count"] = str(count)
    row["telocal_locus"] = locus_key or "."
    row["telocal_active"] = "yes" if count > 0 else "no"

    # Strategy B: refine te_id when the selected one has no expression
    # but another overlapping TE does.
    if count == 0 and hits:
        donor_te_ids = parse_hits_csv(row.get("donor_hits", ""))
        acceptor_te_ids = parse_hits_csv(row.get("acceptor_hits", ""))
        all_te_ids = set(donor_te_ids + acceptor_te_ids)
        for h in hits:
            candidate_key, candidate_count = h[2], h[5]
            if candidate_count > 0:
                # Compare INSERTION ids. This used to take
                # classify_telocal()[0], which is the family ("L1"), and
                # test it against a set of te_ids ("L1PA2_dup1") -- the two
                # vocabularies never intersect, so the refinement could
                # never fire on any key format.
                cand_te_id = telocal_te_id(candidate_key)
                if cand_te_id and cand_te_id in all_te_ids:
                    row["te_id"] = cand_te_id
                    row["telocal_count"] = str(candidate_count)
                    row["telocal_locus"] = candidate_key
                    row["telocal_active"] = "yes"
                    row["te_refined_by_telocal"] = "yes"
                    break

    return row


def main():
    ap = argparse.ArgumentParser(
        description="Annotate chimera junctions with TElocal locus expression."
    )
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--telocal-index", required=True)
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--te-out", required=False,
        help="Optional second output: the gene<->TE subset (direction "
        "gene_to_te / te_to_gene) of --out, same columns.",
    )
    args = ap.parse_args()

    telocal_index = TelocalIndex.load(args.telocal_index)

    with open_read(args.junctions) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        new_cols = [
            "telocal_count", "telocal_locus", "telocal_active",
            "te_refined_by_telocal",
        ]
        out_header = header + new_cols
        rows = []
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            row = annotate_event(row, telocal_index, args.breakpoint_tolerance)
            rows.append(row)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(out_header) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, ".") for c in out_header) + "\n")

    if args.te_out:
        # Same gene<->TE filter classify_chimera_junctions.py's --te-out uses,
        # so the subset is defined identically on both sides.
        te_rows = [
            r for r in rows
            if r.get("direction") in ("gene_to_te", "te_to_gene")
        ]
        os.makedirs(os.path.dirname(args.te_out), exist_ok=True)
        with open_write(args.te_out) as fh:
            fh.write("\t".join(out_header) + "\n")
            for row in te_rows:
                fh.write("\t".join(row.get(c, ".") for c in out_header) + "\n")
        print(f"{args.sample_name}: {len(te_rows)} gene-TE chimeras -> "
              f"{args.te_out}")

    n_active = sum(1 for r in rows if r["telocal_active"] == "yes")
    n_refined = sum(1 for r in rows if r["te_refined_by_telocal"] == "yes")
    print(f"{args.sample_name}: {len(rows)} events, "
          f"{n_active} telocal-active, {n_refined} refined by telocal")


if __name__ == "__main__":
    main()
