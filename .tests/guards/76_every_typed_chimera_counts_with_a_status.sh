#!/usr/bin/env bash
# Guard 76: every typed gene-TE chimera counts, and chimera_status says
# which kind
#
# Only te_initiated / te_terminated / te_exonized used to count, so a pair
# found only as a known TE-driven transcript (annotated_promoter_embedded_te,
# annotated_terminal_exon_embedded_te, annotated_splice) or as antisense
# transcription (antisense_to_gene) dropped out of candidates.tsv.gz. On a
# real 84-sample run that hid e.g. a TE-in-3'-UTR pair a DE analysis on the
# assembly matrix had ranked first. Every typed call now counts toward
# found_by / n_screens, and chimera_status records the kinds, "+"-joined in
# the order novel, annotated, antisense. An untyped event still counts for
# nothing.
#
# Run on its own:   .tests/guards/76_every_typed_chimera_counts_with_a_status.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

python3 - "$T" <<'PY'
import gzip, sys
T = sys.argv[1]
def write(path, cols, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
write(f"{T}/cr.tsv.gz",
      ["event_id", "gene_id", "te_id", "te_subfamily", "te_family", "te_class",
       "canonical", "chimera_type", "n_samples", "total_reads"],
      [["c1", "G_UNTYPED", "TE_U", "L1Md", "L1", "LINE", "no", ".", "1", "3"]])
write(f"{T}/asm.tsv.gz",
      ["transcript_id", "te_id", "te_subfamily", "te_family", "te_class",
       "matched_gene_id", "strand_match", "chimera_type"],
      [["T1", "TE_P", "MTB", "ERVL-MaLR", "LTR", "G_PROM", "yes",
        "annotated_promoter_embedded_te"],
       ["T2", "TE_T", "ID4", "ID", "SINE", "G_TERM", "yes",
        "annotated_terminal_exon_embedded_te"],
       ["T3", "TE_A", "B1_Mm", "Alu", "SINE", "G_ANTI", "no",
        "antisense_to_gene"],
       ["T4", "TE_M", "MT2A", "ERVL", "LTR", "G_MIX", "yes",
        "annotated_promoter_embedded_te"]])
write(f"{T}/sj.tsv.gz",
      ["event_id", "gene_id", "te_id", "te_subfamily", "te_family", "te_class",
       "canonical", "chimera_type", "n_samples", "total_reads"],
      [["s1", "G_SPLICE", "TE_S", "L2", "L2", "LINE", "yes",
        "annotated_splice", "3", "40"],
       ["s2", "G_MIX", "TE_M", "MT2A", "ERVL", "LTR", "yes",
        "te_initiated", "2", "6"]])
PY

if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --assembly "$T/asm.tsv.gz" --sj "$T/sj.tsv.gz" \
      --out "$T/cand.tsv.gz" > "$T/ev.log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev.log"; exit 1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/cand.tsv.gz" --out "$T/table_mqc.json" > "$T/tbl.log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed"; cat "$T/tbl.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, json, sys
T = sys.argv[1]
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False

with gzip.open(f"{T}/cand.tsv.gz", "rt") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    cols = rdr.fieldnames
    cand = {(r["gene_id"], r["te_id"]): r for r in rdr}
check("chimera_status" in cols
      and cols.index("chimera_status") == cols.index("n_screens") + 1,
      f"chimera_status must follow n_screens in candidates.tsv.gz; got {cols}")

want = {
    ("G_PROM", "TE_P"): ("1", "assembly", "annotated"),
    ("G_TERM", "TE_T"): ("1", "assembly", "annotated"),
    ("G_SPLICE", "TE_S"): ("1", "sj", "annotated"),
    ("G_ANTI", "TE_A"): ("1", "assembly", "antisense"),
    ("G_MIX", "TE_M"): ("2", "sj+assembly", "novel+annotated"),
}
for pair, (n, found, status) in want.items():
    r = cand.get(pair)
    check(r is not None, f"{pair} has a typed chimera call and must be a candidate")
    if r is None:
        continue
    got = (r.get("n_screens"), r.get("found_by"), r.get("chimera_status"))
    # found_by order is the script's own; compare as a set of screens
    check(got[0] == n and set(got[1].split("+")) == set(found.split("+"))
          and got[2] == status,
          f"{pair}: want n_screens={n}, found_by={found}, chimera_status={status}; "
          f"got {got}")
check(("G_UNTYPED", "TE_U") not in cand,
      "an untyped chimeric-read event is not a chimera call and must not make a candidate")
log = open(f"{T}/ev.log").read()
check("1 gene-TE pair(s) left out" in log,
      f"the untyped pair must be logged as left out; log: {log!r}")

d = json.load(open(f"{T}/table_mqc.json"))
headers = d.get("headers", {})
check("Status" in headers, f"Candidates table is missing the 'Status' column; got {list(headers)}")
if "Status" in headers and "Types" in headers:
    check(headers["Types"].get("placement") < headers["Status"].get("placement"),
          "Status must sit right after Types")
rows = {k: v for k, v in d.get("data", {}).items()}
# the table shows the top-Screens group: only the 2-screen pair here
mix = rows.get("G_MIX | TE_M", {})
check(mix.get("Status") == "novel+annotated",
      f"table Status for G_MIX | TE_M must be 'novel+annotated'; got {mix.get('Status')!r}")
check(mix.get("Types") == "annotated_promoter_embedded_te,te_initiated",
      f"Types lists every chimera call; got {mix.get('Types')!r}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
