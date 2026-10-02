#!/usr/bin/env bash
# Guard 62: chimera_splice_junctions_igv_bed matches the other two screens
#
# Before this rule existed, chimera.splice_junctions was the only one of the
# three chimera screens with no IGV track at all -- chimera.chimeric_reads and
# chimera.assembly both write one, config-gated by their own
# outputs.write_igv_bed. chimera_splice_junctions_to_igv_bed.py mirrors
# chimera_chimeric_reads_to_igv_bed.py's shape exactly (BED6, same
# gene_to_te/te_to_gene direction filter, same 0-based-half-open coordinate
# math) so the two screens' tracks load and read the same way side by side.
#
# Run on its own:   .tests/guards/62_chimera_splice_junctions_igv_bed_matches_the_other_screens.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- part 1: script-level output ----------------------------------------
# The script runs under Snakemake's `script:` directive (module-level code
# reading snakemake.input/output/params, no main() to call directly) --
# same style as its chimera_chimeric_reads sibling, and the same
# SimpleNamespace-exec technique .tests/unit/test_config_used_mqc.py uses
# for the analogous case.
mkdir -p "$T/igv"
header="event_id\tsample\tchrom\tintron_start\tintron_end\tstrand\tmotif\tcanonical\tannotated\tunique_reads\tmulti_reads\toverhang\tdonor_hits\tacceptor_hits\tdirection\tdirection_ambiguous\tgene_id\tgene_strand\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\tte_initiated_detail\tantisense_flag\tlibrary_strand\ttranscript_strand\tgene_strand_match"
{
  printf "$header\n"
  # kept: gene_to_te, intron [1000, 1050] 1-based inclusive -> BED [999, 1050)
  printf "E1\tS1\tchr1\t1000\t1050\t+\t1\tyes\t1\t7\t0\t20\t.\t.\tgene_to_te\tno\tG1\t+\tT1\t.\t.\t.\tte_exonized\t.\t.\tNA\tNA\tNA\n"
  # kept: te_to_gene, minus strand
  printf "E2\tS1\tchr2\t2000\t2100\t-\t1\tyes\t1\t3\t0\t20\t.\t.\tte_to_gene\tno\tG2\t-\tT2\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\n"
  # dropped: neither gene_to_te nor te_to_gene
  printf "E3\tS1\tchr3\t3000\t3050\t+\t1\tyes\t1\t9\t0\t20\t.\t.\tgene_to_gene\tno\t.\t+\t.\t.\t.\t.\t.\t.\t.\tNA\tNA\tNA\n"
} > "$T/igv/S1_junctions.tsv"

python3 - "$T/igv/S1_junctions.tsv" "$T/igv/S1_junctions.bed" <<'PY' || FAIL=1
import sys, types, os
inp, out = sys.argv[1], sys.argv[2]
script = os.path.join(os.getcwd(), "workflow", "scripts",
                       "chimera_splice_junctions_to_igv_bed.py")
smk = types.SimpleNamespace(input=[inp], output=[out],
                             params=types.SimpleNamespace(sample="S1"))
g = {"snakemake": smk, "__name__": "__main__", "__file__": script}
with open(script) as fh:
    code = fh.read()
exec(compile(code, script, "exec"), g)

lines = open(out).read().splitlines()
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False

check(lines[0].startswith("track name=\"chimera_splice_junctions\""),
      f"track line missing/wrong: {lines[0]!r}")
body = lines[1:]
check(len(body) == 2, f"expected exactly 2 kept rows (E1, E2); got {len(body)}: {body}")
by_event = {row.split("\t")[3]: row.split("\t") for row in body}
check("E3" not in by_event, "E3 (gene_to_gene, neither direction) must be filtered out")
e1 = by_event.get("E1")
check(e1 == ["chr1", "999", "1050", "E1", "7", "+"],
      f"E1 BED row wrong (1-based intron [1000,1050] -> BED [999,1050)): {e1}")
e2 = by_event.get("E2")
check(e2 == ["chr2", "1999", "2100", "E2", "3", "-"],
      f"E2 BED row wrong (1-based intron [2000,2100] -> BED [1999,2100)): {e2}")
sys.exit(0 if ok else 1)
PY

# --- part 2: rule appears/disappears with the outputs.write_igv_bed gate -
# config/test.yaml already has splice_junctions.enabled: true and
# outputs.write_igv_bed: true (set alongside this rule's addition).
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/on.log" 2>&1; then
  echo "ERROR: dry run with write_igv_bed on failed"; tail -30 "$T/on.log"; FAIL=1
elif ! grep -qE "^chimera_splice_junctions_igv_bed\b" "$T/on.log"; then
  echo "ERROR: chimera_splice_junctions_igv_bed not scheduled with write_igv_bed: true"
  FAIL=1
fi

# Anchored on the FULL line (^...$), not a substring match -- a naive range
# ending at the first "write_igv_bed: true" elsewhere (chimeric_reads/
# assembly, earlier in the file) would flip the wrong screen's flag.
sed '/^  splice_junctions:/,/^      write_igv_bed: true$/ s/^      write_igv_bed: true$/      write_igv_bed: false/' \
    config/test.yaml > "$T/off.yaml"
if ! diff -q config/test.yaml "$T/off.yaml" >/dev/null 2>&1; then
  :  # expected: exactly one line changed
else
  echo "ERROR: guard 62 sed did not change anything"; FAIL=1
fi
if ! snakemake --configfile "$T/off.yaml" -n --cores 1 > "$T/off.log" 2>&1; then
  echo "ERROR: dry run with write_igv_bed off failed"; tail -30 "$T/off.log"; FAIL=1
elif grep -qE "^chimera_splice_junctions_igv_bed\b" "$T/off.log"; then
  echo "ERROR: chimera_splice_junctions_igv_bed still scheduled with write_igv_bed: false"
  FAIL=1
fi

exit $FAIL
