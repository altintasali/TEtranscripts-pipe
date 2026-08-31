#!/usr/bin/env bash
# Guard 42: gene_id -> gene_name lookup
#
# Run on its own:   .tests/guards/42_gene_id_gene_name_lookup.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- Every output table is gene_id-keyed, so with an Ensembl or
# GENCODE GTF none of them carry a readable symbol. The lookup is
# how a symbol gets onto any of them; it must be complete (one row
# per gene) and must never leave a blank, or a join silently drops
# or unlabels rows.
mkdir -p "$T/gn"
printf '#!genome-build GRCm38\n' > "$T/gn/ens.gtf"
printf '1\thavana\tgene\t3073253\t3074322\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_name "4933401J01Rik";\n' >> "$T/gn/ens.gtf"
printf '1\thavana\texon\t3073253\t3074322\t.\t+\t.\tgene_id "ENSMUSG00000102693"; transcript_id "T1"; gene_name "4933401J01Rik";\n' >> "$T/gn/ens.gtf"
printf '1\tensembl\tgene\t4807788\t4848410\t.\t+\t.\tgene_id "ENSMUSG00000025900"; gene_name "Rp1";\n' >> "$T/gn/ens.gtf"
# UCSC style: gene_id is already the symbol, no gene_name anywhere
printf 'chr1\tucsc\texon\t100\t200\t.\t+\t.\tgene_id "Actb"; transcript_id "Actb.1";\n' > "$T/gn/ucsc.gtf"
python3 - "$T/gn" <<'PY' || FAIL=1
import builtins, runpy, sys, types
class NS(dict):
    def __getattr__(self, k): return self[k]
T = sys.argv[1]
def run(gtf, out, log):
    builtins.snakemake = types.SimpleNamespace(
        input=NS(gtf=f"{T}/{gtf}"), output=[f"{T}/{out}"], log=[f"{T}/{log}"])
    runpy.run_path("workflow/scripts/gene_name_lookup.py", run_name="__main__")
    import gzip
    rows = [l.rstrip("\n").split("\t") for l in gzip.open(f"{T}/{out}", "rt")]
    return rows[0], dict(rows[1:])
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
header, m = run("ens.gtf", "ens.tsv.gz", "ens.log")
check(header == ["gene_id", "gene_name"], f"unexpected header {header}")
check(len(m) == 2, f"expected 2 genes, got {len(m)}")
check(m.get("ENSMUSG00000102693") == "4933401J01Rik", "symbol not carried across")
check(m.get("ENSMUSG00000025900") == "Rp1", "gene-feature-only symbol missed")
check(all(v for v in m.values()), "a gene_name came out blank")
header, m = run("ucsc.gtf", "ucsc.tsv.gz", "ucsc.log")
check(m.get("Actb") == "Actb", "no-gene_name GTF must fall back to gene_id")
check("WARNING" in open(f"{T}/ucsc.log").read(),
      "a GTF with no gene_name must warn that the lookup is a no-op")
sys.exit(0 if ok else 1)
PY

exit $FAIL
