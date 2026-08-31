#!/usr/bin/env python3
"""Generate the full, rule-level pipeline flowchart.

Reads the DOT graph emitted by `snakemake --rulegraph` from stdin and prints a
Mermaid flowchart with one subgraph per workflow phase. Rule names are mapped
to short friendly labels; any rule not in the map is rendered under its own
name in an "Other" subgraph, so new rules never break the generator.

This is the exhaustive, engineer-facing diagram (one node per rule) -- it
lives in docs/pipeline-flowchart.md (mirrored by hand into the wiki's
"Full Pipeline Flowchart" page) rather than the README, which instead has a
small, hand-authored, conceptual diagram aimed at a general user.

Usage:
    snakemake --configfile config/test.yaml --rulegraph | workflow/scripts/flowchart.py
    snakemake --configfile config/test.yaml --rulegraph | workflow/scripts/flowchart.py --update-doc

--update-doc rewrites the block between the `<!-- flowchart:start -->` and
`<!-- flowchart:end -->` markers in docs/pipeline-flowchart.md in place. The
CI workflow runs this and fails on a diff, so the committed diagram can't go
stale (the wiki's copy of it is a manual refresh, same as the rest of the
wiki -- not CI-enforced).
"""
import os
import re
import sys

DOC = "docs/pipeline-flowchart.md"
DOC_HEADER = (
    "# Full Pipeline Flowchart\n\n"
    "Auto-generated from `snakemake --rulegraph` by `workflow/scripts/"
    "flowchart.py` -- do not hand-edit (regenerate instead, see the script's "
    "docstring). One node per rule; for a simpler, conceptual diagram see the "
    "main [README](https://github.com/altintasali/TEtranscripts-pipe#readme).\n\n"
)
FLOW_START = "<!-- flowchart:start -->"
FLOW_END = "<!-- flowchart:end -->"

# rule -> (friendly label, phase). Insertion order sets the node order within
# a phase; PHASES sets the subgraph order. Rules missing from this map still
# get a node (their own name, phase "Other").
LABELS = {
    "gunzip_reference": ("gunzip reference (if .gz refs)", "Reference (once)"),
    "star_index": ("STAR index", "Reference (once)"),
    "gtf_to_genepred": ("GTF -> genePred", "Reference (once)"),
    "genepred_to_bed12": ("genePred -> BED12", "Reference (once)"),
    "cat_fastq": ("concat lanes", "Per sample"),
    "trim_galore_pe": ("Trim Galore! (paired)", "Per sample"),
    "trim_galore_se": ("Trim Galore! (single-end)", "Per sample"),
    "star_align": ("STAR align", "Per sample"),
    "samtools_sort": ("samtools sort", "Per sample"),
    "samtools_index": ("samtools index", "Per sample"),
    "fastqc_raw": ("FastQC (raw)", "Per sample"),
    "rseqc_infer_experiment": ("RSeQC infer_experiment", "Per sample"),
    "determine_strandedness": ("determine strandedness", "Per sample"),
    "tecount": ("TEcount", "Quantification + QC"),
    "tetranscripts_diffexp": ("TEtranscripts + DESeq2", "Quantification + QC"),
    "tecount_counts": ("tecount counts matrix", "Quantification + QC"),
    "tecount_qc_transform": ("sample-QC transform (vst/rlog/log2)", "Quantification + QC"),
    "tecount_qc": ("sample-QC plots (PCA + clustering)", "Quantification + QC"),
    "tecount_summary": ("tecount summary barplots (assignment + TE class)", "Quantification + QC"),
    "telocal_locind": ("TElocal locus index", "Reference (once)"),
    "telocal": ("TElocal", "Quantification + QC"),
    "telocal_counts": ("telocal counts matrix", "Quantification + QC"),
    "telocal_qc_transform": ("telocal sample-QC transform (log2/vst/rlog)", "Quantification + QC"),
    "telocal_qc": ("telocal sample-QC plots (PCA + clustering)", "Quantification + QC"),
    "telocal_summary": ("telocal summary barplots (assignment + TE class)", "Quantification + QC"),
    "cleanup_star_index": ("remove STAR index (if keep: false)", "Reference (once)"),
    "cleanup_telocal_index": ("remove TElocal index (if keep: false)", "Quantification + QC"),
    "gene_name_lookup": ("gene_id -> gene_name lookup", "Reference (once)"),
    "annotation_to_bed": ("annotation -> BED tracks", "Chimera screen"),
    "parse_chimeric_junctions": ("parse chimeric junctions", "Chimera screen"),
    "chimera_telocal_annotate": ("annotate junctions with TElocal counts", "Chimera screen"),
    "junction_qc": ("junction QC", "Chimera screen"),
    "junction_qc_barplot": ("junction QC barplot", "Chimera screen"),
    "junction_highlights": ("read-screen notes (blind spot + counts)", "Chimera screen"),
    "chimera_evidence": ("unified gene-TE evidence catalogue", "Chimera screen"),
    "chimera_evidence_heatmap": ("evidence correlation + candidate heatmaps", "Chimera screen"),
    "chimera_evidence_guide": ("how to weigh the evidence + composition", "Chimera screen"),
    "sample_evidence_status": ("per-sample evidence status grid", "Chimera screen"),
    "chimera_igv_bed": ("IGV BED track", "Chimera screen"),
    "chimera_counts": ("chimera counts matrix", "Chimera screen"),
    "sample_qc_transform": ("sample-QC transform", "Chimera screen"),
    "sample_qc": ("sample-QC plots", "Chimera screen"),
    "software_versions": ("software versions", "Quantification + QC"),
    "config_used": ("config used", "Quantification + QC"),
    "evidence_overview": ("evidence overview ('start here')", "Quantification + QC"),
    "strandedness_check": ("strandedness check (declared vs inferred)", "Quantification + QC"),
    "benchmark_summary": ("resource-usage summary", "Quantification + QC"),
    "multiqc": ("MultiQC", "Quantification + QC"),
}
PHASES = ["Reference (once)", "Per sample", "Quantification + QC", "Chimera screen", "Other"]
IGNORED_RULES = {"all"}
# benchmark_summary depends on every rule's benchmark file, and multiqc on
# most rules' outputs -- drawing all those fan-in edges would turn the
# diagram into a spiderweb. Keep the nodes but only representative edges
# (target -> allowed sources; empty tuple suppresses all fan-in).
AGGREGATOR_FAN_IN = {
    "benchmark_summary": (),
    "multiqc": ("software_versions", "benchmark_summary"),
}

NODE_RE = re.compile(r'^\s*(\d+)\[label\s*=\s*"([^"]+)"')
EDGE_RE = re.compile(r"^\s*(\d+)\s*->\s*(\d+)")


def parse(dot):
    names = {}
    edges = []
    for line in dot.splitlines():
        m = NODE_RE.match(line)
        if m:
            names[m.group(1)] = m.group(2)
            continue
        m = EDGE_RE.match(line)
        if m:
            edges.append((m.group(1), m.group(2)))
    return names, edges


def mermaid_block(names, edges):
    names = {rid: r for rid, r in names.items() if r not in IGNORED_RULES}
    rule_ids = {r: rid for rid, r in names.items()}
    used = set(names.values())

    phase_order = {p: [] for p in PHASES}
    for rule in LABELS:
        if rule in used:
            label, phase = LABELS[rule]
            phase_order[phase].append((rule, label))
    for rule in sorted(used - set(LABELS)):
        phase_order["Other"].append((rule, rule))

    lines = ["```mermaid", "flowchart LR"]
    slug = {p: re.sub(r"[^a-z]+", "_", p.lower()).strip("_") for p in PHASES}
    for phase in PHASES:
        if not phase_order[phase]:
            continue
        lines.append(f'    subgraph {slug[phase]}["{phase}"]')
        for rule, label in phase_order[phase]:
            lines.append(f'        {rule}["{label}"]')
        lines.append("    end")
    for a, b in sorted(
        (names[a], names[b])
        for a, b in edges
        if a in names and b in names and b not in IGNORED_RULES
    ):
        allowed = AGGREGATOR_FAN_IN.get(b)
        if allowed is None or a in allowed:
            lines.append(f"    {a} --> {b}")
    lines.append("```")
    return "\n".join(lines)


def update_doc(block):
    if os.path.isfile(DOC):
        with open(DOC) as fh:
            text = fh.read()
    else:
        text = DOC_HEADER + FLOW_START + "\n" + FLOW_END + "\n"
        os.makedirs(os.path.dirname(DOC), exist_ok=True)

    start = text.find(FLOW_START)
    end = text.find(FLOW_END)
    if start == -1 or end == -1 or end <= start:
        sys.exit(
            f"error: {DOC} is missing the {FLOW_START!r} / {FLOW_END!r} "
            "markers around the flowchart"
        )
    end = end + len(FLOW_END)
    with open(DOC, "w") as fh:
        fh.write(text[:start] + FLOW_START + "\n" + block + "\n" + FLOW_END + text[end:])


def main():
    block = mermaid_block(*parse(sys.stdin.read()))
    if "--update-doc" in sys.argv:
        update_doc(block)
    else:
        print(block)


if __name__ == "__main__":
    main()
