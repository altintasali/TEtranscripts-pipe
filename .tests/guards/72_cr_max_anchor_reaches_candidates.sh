#!/usr/bin/env bash
# Guard 72: chimeric-read anchor length reaches the candidates
#
# The chimeric-read counterpart of SJ.out.tab's overhang: per read, the
# aligned length of its SHORTER segment (M bases in Chimeric.out.junction's
# CIGAR columns 12 / 14). On a real run local gene-TE events had a median
# best anchor of ~26 bp vs ~18 bp for trans / far ones.
#   - classify_chimera_chimeric_reads.py: max_anchor = best read per event;
#   - chimera_chimeric_reads_counts.py: max across samples ("." when no
#     table carries the column);
#   - chimera_evidence.py: cr_max_anchor from the pair's CR CALLS only;
#   - the Candidates table: "CR max anchor" in the CR block.
#
# Run on its own:   .tests/guards/72_cr_max_anchor_reaches_candidates.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# GENE1 (+) exon [1000,1200), TE_UP [500,700): a TE->gene junction, as in
# guard 67 (te_initiated on any library).
printf 'chr1\t1000\t1200\tGENE1\t.\t+\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\n' > "$T/exons.bed"
printf 'chr1\t500\t700\tTE_UP\t.\t+\tERVL-MaLR\tLTR\tMTA\n' > "$T/te.bed"
# two reads on one event: shorter segments 25 and 35 -> max_anchor 35
{
  printf 'chr1\t700\t+\tchr1\t1000\t+\t1\t0\t0\tr1\t650\t50M25S\t1001\t50S25M\n'
  printf 'chr1\t700\t+\tchr1\t1000\t+\t1\t0\t0\tr2\t660\t40M35S\t1001\t40S35M\n'
} > "$T/S1.junction"

if ! python3 workflow/scripts/classify_chimera_chimeric_reads.py \
      --junctions "$T/S1.junction" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --sample S1 --library-strandedness no \
      --out "$T/S1.tsv.gz" > "$T/cls.log" 2>&1; then
  echo "ERROR: classify_chimera_chimeric_reads.py failed"; cat "$T/cls.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, subprocess, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
def read(p):
    with gzip.open(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
def write(p, rows, cols):
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

s1 = read(f"{T}/S1.tsv.gz")
check(len(s1) == 1 and s1[0].get("max_anchor") == "35",
      f"classifier: max_anchor must be the best read's shorter segment (35); got "
      f"{[r.get('max_anchor') for r in s1]}")
check(s1[0].get("chimera_type") == "te_initiated", f"fixture must be a call; got {s1[0].get('chimera_type')!r}")
cols = list(s1[0].keys())

# S2: same event, better read (50); S3: an older table without the column
s2 = [dict(s1[0], sample="S2", max_anchor="50")]
write(f"{T}/S2.tsv.gz", s2, cols)
old_cols = [c for c in cols if c != "max_anchor"]
write(f"{T}/S3.tsv.gz", [dict(s1[0], sample="S3")], old_cols)

def counts(tables, names, out):
    subprocess.run([sys.executable, "workflow/scripts/chimera_chimeric_reads_counts.py",
                    "--tables", *tables, "--sample-names", *names,
                    "--out-events", f"{T}/{out}_all.tsv.gz",
                    "--out-te-events", f"{T}/{out}_te.tsv.gz"],
                   check=True, capture_output=True)
    return read(f"{T}/{out}_te.tsv.gz")

m = counts([f"{T}/S1.tsv.gz", f"{T}/S2.tsv.gz"], ["S1", "S2"], "m")
check(m[0].get("max_anchor") == "50",
      f"counts: max_anchor is the max across samples (50), not the first sample's; "
      f"got {m[0].get('max_anchor')!r}")
o = counts([f"{T}/S3.tsv.gz"], ["S3"], "o")
check(o[0].get("max_anchor") == ".",
      f"counts: a table without the column gives '.', not 0; got {o[0].get('max_anchor')!r}")

# evidence: the call (anchor 50) plus a non-call event of the same pair with
# a bigger anchor, which must not count
merged = m + [dict(m[0], event_id="x", chimera_type="antisense_to_gene", max_anchor="99")]
write(f"{T}/cr.tsv.gz", merged, list(m[0].keys()))
subprocess.run([sys.executable, "workflow/scripts/chimera_evidence.py",
                "--junction", f"{T}/cr.tsv.gz", "--out", f"{T}/cand.tsv.gz"],
               check=True, capture_output=True)
cand = read(f"{T}/cand.tsv.gz")
check(len(cand) == 1 and cand[0].get("cr_max_anchor") == "50",
      f"evidence: cr_max_anchor from the CALL only (50, not the non-call's 99); "
      f"got {[r.get('cr_max_anchor') for r in cand]}")

subprocess.run([sys.executable, "workflow/scripts/chimera_candidates_table_mqc.py",
                "--evidence", f"{T}/cand.tsv.gz", "--out", f"{T}/table.json"],
               check=True, capture_output=True)
import json
d = json.load(open(f"{T}/table.json"))
e = d["data"].get("GENE1 | TE_UP", {})
check("CR max anchor" in d.get("headers", {}) and e.get("CR max anchor") == 50,
      f"table: 'CR max anchor' column with 50; got {e}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
