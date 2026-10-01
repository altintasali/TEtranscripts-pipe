#!/usr/bin/env bash
# Guard 71: SJ mapping quality reaches the candidates
#
# Reads from young TE families map to many copies, and spurious junctions
# arise there. STAR's SJ.out.tab already reports unique_reads, multi_reads
# and overhang per junction per sample; none of it reached the candidates:
#   - chimera_splice_junctions_counts.py copied multi_reads / overhang from
#     the FIRST sample that saw a junction; it now sums multi_reads and takes
#     the max overhang across samples (unique reads were already summed, as
#     total_reads);
#   - chimera_evidence.py adds sj_unique_fraction = unique / (unique + multi)
#     and sj_max_overhang, from the pair's SJ chimera CALLS only;
#   - the Candidates table shows both in its SJ block.
#
# Run on its own:   .tests/guards/71_sj_mapping_quality_reaches_candidates.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

hdr="event_id\tsample\tchrom\tintron_start\tintron_end\tstrand\tmotif\tcanonical\tannotated\tunique_reads\tmulti_reads\toverhang\tdonor_hits\tacceptor_hits\tdirection\tdirection_ambiguous\tgene_id\tgene_strand\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\tte_initiated_detail\tantisense_flag\tlibrary_strand\ttranscript_strand\tgene_strand_match"
row() {  # sample unique multi overhang
  printf "chr1:701:1000:+\t%s\tchr1\t701\t1000\t+\t1\tyes\t0\t%s\t%s\t%s\t.\t.\tte_to_gene\tno\tG1\t+\tTE1\tsf\tfam\tLTR\tte_initiated\tupstream\t.\tno\t+\tyes\n" "$1" "$2" "$3" "$4"
}
{ printf "%b\n" "$hdr"; row S1 3 1 20; } | gzip -c > "$T/S1.tsv.gz"
{ printf "%b\n" "$hdr"; row S2 5 3 35; } | gzip -c > "$T/S2.tsv.gz"

if ! python3 workflow/scripts/chimera_splice_junctions_counts.py \
      --tables "$T/S1.tsv.gz" "$T/S2.tsv.gz" --sample-names S1 S2 \
      --out-events "$T/all.tsv.gz" --out-te-events "$T/te.tsv.gz" > "$T/counts.log" 2>&1; then
  echo "ERROR: chimera_splice_junctions_counts.py failed"; cat "$T/counts.log"; exit 1
fi

# merged SJ table for chimera_evidence: G1/TE1 has the real call (8 unique,
# 4 multi, overhang 35 -- the counts step's own output, appended below) plus
# a non-call annotated splice with heavy multi-mapping that must NOT enter
# the fraction; G2/TE2 has no SJ call at all.
python3 - "$T" <<'PY'
import csv, gzip, sys
T = sys.argv[1]
with gzip.open(f"{T}/te.tsv.gz", "rt") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
cols = list(rows[0].keys())
extra = dict(rows[0], event_id="annot", chimera_type="annotated_splice",
             total_reads="100", multi_reads="900", overhang="99")
with gzip.open(f"{T}/sj_merged.tsv.gz", "wt") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows + [extra]:
        fh.write("\t".join(r[c] for c in cols) + "\n")
with gzip.open(f"{T}/cr.tsv.gz", "wt") as fh:
    fh.write("event_id\tgene_id\tte_id\tchimera_type\tn_samples\ttotal_reads\n")
    fh.write("c1\tG2\tTE2\tte_initiated\t1\t2\n")
PY

if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --sj "$T/sj_merged.tsv.gz" --out "$T/cand.tsv.gz" > "$T/ev.log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev.log"; exit 1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/cand.tsv.gz" --out "$T/table.json" > "$T/tbl.log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed"; cat "$T/tbl.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, json, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
def rows(p):
    with gzip.open(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

ev = rows(f"{T}/te.tsv.gz")[0]
check((ev["total_reads"], ev["multi_reads"], ev["overhang"]) == ("8", "4", "35"),
      f"cohort SJ event must sum unique/multi reads and take the max overhang "
      f"across samples (want 8/4/35); got {ev['total_reads']}/{ev['multi_reads']}/{ev['overhang']}")

cand = {(r["gene_id"], r["te_id"]): r for r in rows(f"{T}/cand.tsv.gz")}
g1 = cand.get(("G1", "TE1"), {})
check(g1.get("sj_unique_fraction") == "0.667" and g1.get("sj_max_overhang") == "35",
      f"sj_unique_fraction = 8/(8+4) and sj_max_overhang = 35 from the CALL only "
      f"(the annotated splice's 900 multi reads / overhang 99 must be ignored); got "
      f"{g1.get('sj_unique_fraction')!r}/{g1.get('sj_max_overhang')!r}")
g2 = cand.get(("G2", "TE2"), {})
check(g2.get("sj_unique_fraction") == "." and g2.get("sj_max_overhang") == ".",
      f"no SJ call -> both '.'; got {g2.get('sj_unique_fraction')!r}/{g2.get('sj_max_overhang')!r}")

d = json.load(open(f"{T}/table.json"))
e = d["data"].get("G1 | TE1", {})
check(e.get("SJ unique fraction") == 0.667 and e.get("SJ max overhang") == 35,
      f"Candidates table must show both in its SJ block; got {e}")
e2 = d["data"].get("G2 | TE2", {})
check("SJ unique fraction" not in e2 and "SJ max overhang" not in e2,
      f"no SJ call -> blank cells, not 0; got {e2}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
