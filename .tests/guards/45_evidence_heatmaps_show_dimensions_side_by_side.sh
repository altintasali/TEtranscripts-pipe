#!/usr/bin/env bash
# Guard 45: evidence heatmaps show dimensions side by side
#
# Run on its own:   .tests/guards/45_evidence_heatmaps_show_dimensions_side_by_side.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The whole point of these plots is that they expose a
# dimension pointing the WRONG way, which a combined score hides.
# Fixture bakes in the anti-correlation measured on real data (TE
# locus expression vs splice motif); the correlation heatmap must
# recover it as a clearly negative cell.
mkdir -p "$T/hm/qc"
python3 - "$T/hm" <<'PY'
import gzip, random, sys
T = sys.argv[1]; random.seed(7)
cols = ["gene_id","te_id","te_subfamily","te_family","te_class","found_by",
      "evidence","n_evidence","junction_events","junction_reads","junction_max_samples",
      "junction_canonical","junction_chimera_types","telocal_active",
      "assembly_transcripts","assembly_chimera_types","assembly_strand_match",
      "assembly_transcript_ids"]
with gzip.open(f"{T}/ev.tsv.gz","wt") as fh:
  fh.write("\t".join(cols)+"\n")
  for i in range(3000):
      canon = random.random() < 0.086
      active = "yes" if random.random() < (0.80 if canon else 0.93) else "no"
      fb = "both" if random.random() < 0.0042 else ("junction" if random.random() < 0.6 else "assembly")
      asm = random.randint(1,14) if fb in ("both","assembly") else 0
      jr = random.randint(1,40) if fb in ("both","junction") else 0
      nsamp = (random.choice([1,1,1,1,2,3]) if jr else 0)
      strand = "yes" if asm and random.random()<0.7 else "."
      flags = ([f for f, on in (("canonical", canon),
                                ("multi_sample", nsamp > 1),
                                ("both_screens", fb == "both"),
                                ("assembly_strand_match", strand == "yes"))
                if on])
      r = {"gene_id":f"ENSMUSG{i:011d}","te_id":f"L1x_dup{i}","te_subfamily":"L1x",
           "te_family":"L1","te_class":"LINE","found_by":fb,
           # mirrors chimera_evidence.py's flag set so the fixture
           # stays the same shape as the real table
           "evidence":",".join(flags) or ".","n_evidence":len(flags),
           "junction_events":random.randint(1,3) if jr else 0,"junction_reads":jr,
           "junction_max_samples":nsamp,
           "junction_canonical":"yes" if canon else "no",
           "junction_chimera_types":"te_exonized","telocal_active":active if jr else ".",
           "assembly_transcripts":asm,"assembly_chimera_types":"te_exonized",
           "assembly_strand_match":strand,
           "assembly_transcript_ids":"MSTRG.1"}
      fh.write("\t".join(str(r[c]) for c in cols)+"\n")
with gzip.open(f"{T}/sym.tsv.gz","wt") as fh:
  fh.write("gene_id\tgene_name\n")
  for i in range(3000): fh.write(f"ENSMUSG{i:011d}\tGene{i}\n")
PY
if ! python3 workflow/scripts/chimera_evidence_heatmap.py \
  --evidence "$T/hm/ev.tsv.gz" --gene-names "$T/hm/sym.tsv.gz" \
  --out-correlation "$T/hm/qc/chimera_evidence_correlation_mqc.json" \
  --out-candidates "$T/hm/qc/chimera_evidence_candidates_mqc.json" \
  > "$T/hm/log" 2>&1; then
echo "ERROR: chimera_evidence_heatmap.py failed"; cat "$T/hm/log"; FAIL=1
else
python3 - "$T/hm" <<'PY' || FAIL=1
import json, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
d = json.load(open(f"{T}/qc/chimera_evidence_correlation_mqc.json"))
check(d["plot_type"] == "heatmap", "correlation must be a heatmap")
check("data" in d and "xcats" in d and "ycats" in d, "heatmap needs data/xcats/ycats")
labs = d["xcats"]
n = len(labs)
check(len(d["data"]) == n and all(len(r) == n for r in d["data"]), "matrix not square")
check(all(abs(d["data"][i][i] - 1.0) < 1e-6 for i in range(n)), "diagonal must be 1.0")
i, j = labs.index("Splice motif"), labs.index("TE locus expressed")
check(abs(d["data"][i][j] - d["data"][j][i]) < 1e-9, "matrix must be symmetric")
check(d["data"][i][j] < -0.05,
      f"anti-correlation baked into the fixture not recovered: {d['data'][i][j]}")
c = json.load(open(f"{T}/qc/chimera_evidence_candidates_mqc.json"))
check(c["plot_type"] == "heatmap", "candidates must be a heatmap")
check(len(c["ycats"]) == len(c["data"]), "row labels and rows disagree")
check(all(len(r) == n for r in c["data"]), "candidate rows must cover every dimension")
check(all(0 <= v <= 100 for r in c["data"] for v in r), "cells must be percentiles")
check(any("/" in y for y in c["ycats"]), "rows must be labelled gene / te_id")
# MultiQC labels every heatmap row a "sample" in its own toolbox, so the axes
# must say what they actually are -- 7 evidence types read as 7 samples.
check(d["pconfig"].get("xlab") == "Evidence type",
      f"correlation x-axis must name evidence types, got {d['pconfig'].get('xlab')!r}")
check(d["pconfig"].get("ylab") == "Evidence type",
      f"correlation y-axis must name evidence types, got {d['pconfig'].get('ylab')!r}")
check(c["pconfig"].get("ylab") == "Gene / TE pair",
      f"candidate rows are pairs, not samples; got {c['pconfig'].get('ylab')!r}")
check("evidence types" in d.get("description", "").lower(),
      "the correlation section must say its axes are evidence types, not samples")
sys.exit(0 if ok else 1)
PY
fi
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
  -o "$T/hm/out" -n r "$T/hm/qc" > "$T/hm/render.log" 2>&1; then
echo "ERROR: multiqc failed on the evidence heatmaps"
tail -30 "$T/hm/render.log"; FAIL=1
elif ! grep -q "Evidence structure - correlation" "$T/hm/out/r.html"; then
echo "ERROR: evidence correlation heatmap not rendered"; FAIL=1
elif ! grep -q "Evidence structure - leaders by dimension" "$T/hm/out/r.html"; then
echo "ERROR: candidate heatmap not rendered"; FAIL=1
fi

exit $FAIL
