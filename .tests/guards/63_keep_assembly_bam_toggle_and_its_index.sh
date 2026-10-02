#!/usr/bin/env bash
# Guard 63: outputs.keep_assembly_bam gates chimera_assembly_bam_index, and
# both the BAM and its index are temp()-wrapped together when it is false.
#
# chimera.assembly's own STAR pass (star_align_for_assembly) is a SEPARATE,
# private alignment -- the main results/star/ BAM everything else uses is
# untouched (see that rule's own comment for why). Nothing reads this
# private BAM after stringtie_assemble/stringtie_requantify consume it, so
# it is a real, avoidable disk cost when kept. outputs.keep_assembly_bam
# (default true) controls whether it -- and its samtools index, added
# alongside it -- survive the run; both use the same _maybe_temp() helper
# outputs.keep_merged_fastq/keep_trimmed_fastq already use.
#
# Run on its own:   .tests/guards/63_keep_assembly_bam_toggle_and_its_index.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- part 1: source-level wiring -- both outputs go through _maybe_temp() ---
if ! grep -q 'aln=_maybe_temp(' workflow/rules/chimera_assembly.smk; then
  echo "ERROR: star_align_for_assembly's aln output is not wrapped in _maybe_temp()"
  FAIL=1
fi
if ! sed -n '/^rule chimera_assembly_bam_index:/,/^rule /p' workflow/rules/chimera_assembly.smk \
     | grep -q '_maybe_temp('; then
  echo "ERROR: chimera_assembly_bam_index's output is not wrapped in _maybe_temp()"
  FAIL=1
fi

# --- part 2: schema declares the toggle -------------------------------
python3 - <<'PY' || FAIL=1
import yaml, sys
schema = yaml.safe_load(open("workflow/schemas/config.schema.yaml"))
prop = schema["properties"]["outputs"]["properties"].get("keep_assembly_bam")
if prop is None or prop.get("type") != "boolean":
    print("ERROR: outputs.keep_assembly_bam missing or not boolean in config.schema.yaml")
    sys.exit(1)
PY

# --- part 3: default (true) schedules the index; false does not --------
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/on.log" 2>&1; then
  echo "ERROR: dry run with keep_assembly_bam at its default failed"; tail -30 "$T/on.log"; FAIL=1
elif ! grep -qE "^chimera_assembly_bam_index\b" "$T/on.log"; then
  echo "ERROR: chimera_assembly_bam_index not scheduled with keep_assembly_bam defaulting true"
  FAIL=1
fi

python3 - "$T/off.yaml" <<'PY'
import sys
content = open("config/test.yaml").read()
idx = content.index("outputs:")
end = content.index("\n\n", idx)
block = content[idx:end]
open(sys.argv[1], "w").write(content[:idx] + block + "\n  keep_assembly_bam: false" + content[end:])
PY
if ! grep -q "keep_assembly_bam: false" "$T/off.yaml"; then
  echo "ERROR: guard 63 did not write the keep_assembly_bam: false override"; FAIL=1
elif ! snakemake --configfile "$T/off.yaml" -n --cores 1 > "$T/off.log" 2>&1; then
  echo "ERROR: dry run with keep_assembly_bam: false failed"; tail -30 "$T/off.log"; FAIL=1
elif grep -qE "^chimera_assembly_bam_index\b" "$T/off.log"; then
  echo "ERROR: chimera_assembly_bam_index still scheduled with keep_assembly_bam: false -- indexing a BAM about to be temp()-deleted"
  FAIL=1
fi
# star_align_for_assembly itself must still run either way -- only the
# KEEPING of its output (and whether it gets indexed) is gated, not the
# alignment -- stringtie_assemble/stringtie_requantify need it regardless.
for log in "$T/on.log" "$T/off.log"; do
  if ! grep -qE "^star_align_for_assembly\b" "$log"; then
    echo "ERROR: star_align_for_assembly not scheduled in $log -- the BAM must still be produced for StringTie either way"
    FAIL=1
  fi
done

# --- part 4: config_used_mqc.py renders the resolved value -------------
python3 - "$T" <<'PY' || FAIL=1
import sys, types, os, json
T = sys.argv[1]
sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
import config_used_mqc

base_config = {
    "samples": {}, "ref": {"gtf": "g.gtf", "te_gtf": "te.gtf"},
    "star": {}, "trimming": {}, "strandedness": {"min_fraction": 0.8},
    "tetranscripts": {"mode": "multi"}, "chimera": {},
}
base_params = {
    "_samples": ["S1"], "_sample_count": 1, "_sjdb_overhang": "auto",
    "_star_index": "idx", "_trim_enabled": True,
    "_tecount_qc_enabled": True, "_tecount_qc": {},
    "_chimera_chimeric_reads_enabled": True,
    "_telocal_enabled": True, "_telocal_locind_auto": True,
    "_telocal_qc_enabled": True, "_telocal_qc": {},
    "_keep_merged_fastq": True, "_keep_trimmed_fastq": True,
    "_keep_star_index": True, "_keep_telocal_index": True,
}
ok = True
for keep in (True, False):
    params = dict(base_params, _keep_assembly_bam=keep)
    out_path = os.path.join(T, f"config_used_{keep}.json")
    smk = types.SimpleNamespace(config=base_config, params=params,
                                 log=os.path.join(T, f"config_used_{keep}.log"),
                                 output=out_path)
    config_used_mqc.main(smk)
    doc = json.load(open(out_path))
    row = doc["data"].get(config_used_mqc.row_label("outputs.keep_assembly_bam"), {})
    got = row.get("value")
    if got != str(keep):
        print(f"ERROR: outputs.keep_assembly_bam rendered {got!r}, expected {str(keep)!r}")
        ok = False
sys.exit(0 if ok else 1)
PY

exit $FAIL
