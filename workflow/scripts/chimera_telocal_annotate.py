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
      --telocal-tables {sample}.cntTable.gz ... \\
      --sample-name {sample} \\
      --breakpoint-tolerance 0 \\
      --out {sample}_junctions_telocal.tsv.gz
"""
import argparse
import bisect
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def parse_telocal_locus(key):
    """Extract (chrom, start, end) from a TElocal locus key.

    TElocal TE keys have the form
      chrom:start:end(family:strand):gene_id:family_id:class_id
    Gene keys (no parentheses in transcript portion) return None.
    """
    m = re.match(r"^(.+):(\d+):(\d+)\(", key)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def classify_telocal(key):
    """Return (family, class) from a TElocal TE key, or (None, None) for genes.

    TElocal TE keys have the form:
      chrom:start:end(family:strand):gene_id:family:class
    The last two colon-separated fields after the locus coordinates are
    family and class.  Gene entries (no parentheses) return (None, None).
    """
    m = re.match(r"^(.+):(\d+):(\d+)\(", key)
    if not m:
        return None, None
    locus_end = m.end()
    suffix = key[locus_end:]
    # suffix is like "):L1PA2:L1:LINE/L1" or "):AluYb:SINE"
    # strip leading ')' or '):'
    suffix = suffix.lstrip("):").strip(":")
    parts = suffix.split(":")
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    if len(parts) == 1:
        return parts[0], None
    return None, None


def load_telocal_counts(paths):
    """Build a per-chromosome interval index from TElocal cntTables.

    Returns {chrom: (feats, max_end)}, where `feats` is a list of
    (start, end, locus_key, family, cls, count) sorted by start, and
    `max_end[i]` is the running maximum end coordinate over feats[0..i] --
    used by overlapping_telocal() for a correct (not just "usually correct")
    interval overlap query, mirroring parse_chimeric_junctions.py's
    load_bed()/overlapping(). count is summed across all input files
    (samples).
    """
    raw = {}
    for path in paths:
        with open_read(path) as fh:
            fh.readline()
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                key = parts[0]
                try:
                    count = int(parts[1])
                except ValueError:
                    count = 0
                coords = parse_telocal_locus(key)
                if coords is None:
                    continue
                chrom, start, end = coords
                family, cls = classify_telocal(key)
                bucket = raw.setdefault(chrom, {})
                if key not in bucket:
                    bucket[key] = (start, end, key, family, cls, count)
                else:
                    s, e, k, f, c, n = bucket[key]
                    bucket[key] = (s, e, k, f, c, n + count)

    tracks = {}
    for chrom, d in raw.items():
        feats = sorted(d.values(), key=lambda x: (x[0], x[1]))
        running = float("-inf")
        max_end = []
        for f in feats:
            running = max(running, f[1])
            max_end.append(running)
        tracks[chrom] = (feats, max_end)
    return tracks


def overlapping_telocal(track, chrom, start0, end0):
    """TElocal loci overlapping [start0, end0) on *chrom*.

    Correct for arbitrarily long/nested/overlapping loci, not just the
    common case -- see parse_chimeric_junctions.py's overlapping() for the
    rationale behind bisecting on the running max_end rather than start.
    """
    if chrom not in track:
        return []
    feats, max_end = track[chrom]
    lo = bisect.bisect_right(max_end, start0)
    hits = []
    for i in range(lo, len(feats)):
        s, e, key, family, cls, count = feats[i]
        if s >= end0:
            break
        if e > start0:
            hits.append((s, e, key, family, cls, count))
    return hits


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

    hits = overlapping_telocal(telocal_index, te_chrom, q_start, q_end)
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
                cand_family = classify_telocal(candidate_key)[0]
                if cand_family and cand_family in all_te_ids:
                    row["te_id"] = cand_family
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
    ap.add_argument("--telocal-tables", required=True, nargs="+")
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    telocal_index = load_telocal_counts(args.telocal_tables)

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

    n_active = sum(1 for r in rows if r["telocal_active"] == "yes")
    n_refined = sum(1 for r in rows if r["te_refined_by_telocal"] == "yes")
    print(f"{args.sample_name}: {len(rows)} events, "
          f"{n_active} telocal-active, {n_refined} refined by telocal")


if __name__ == "__main__":
    main()
