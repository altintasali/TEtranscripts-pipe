# -----------------------------------------------------------------------------
# Software versions section for MultiQC.
#
# Writes the pinned tool versions (config["versions"]) in the format
# MultiQC's "Software Versions" section expects (a `*_mqc_versions.yaml` /
# `.yml` file, per MultiQC 1.33's search pattern -- it scans for that
# filename, not `multiqc_software_versions.yml`). MultiQC then shows the
# versions in the Software Versions report section and tags the STAR, FastQC,
# samtools etc. rows in the General Stats table with them, so a report always
# records exactly which releases each step ran with.
# -----------------------------------------------------------------------------
rule software_versions:
    output:
        "results/versions/rnaseq_mqc_versions.yml",
    threads: get_resources("software_versions")["threads"]
    resources:
        mem_mb=get_resources("software_versions")["mem_mb"],
        runtime=get_resources("software_versions")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/software_versions/software_versions.txt",
    log:
        "results/pipeline_info/logs/software_versions.log",
    run:
        import yaml

        # Every pinned tool, grouped by pipeline stage -- MultiQC renders
        # each top-level YAML key as its own header in the Software
        # Versions section, so splitting into stage groups (instead of one
        # flat "TEtranscripts-pipe" list) makes the section scannable as
        # the tool count grows. Tools whose versions MultiQC auto-detects
        # from logs (STAR, FastQC, samtools...) are listed here too, so the
        # section stays complete even when a log omits them.
        versions = {
            "Trimming & QC": {
                "Trim Galore!": V["trim_galore"],
                "FastQC": V["fastqc"],
                "MultiQC": V["multiqc"],
            },
            "Alignment": {
                "STAR": V["star"],
                "samtools": V["samtools"],
                "RSeQC": V["rseqc"],
            },
            "Annotation tools": {
                "UCSC tools": f"{V['ucsc_gtftogenepred']}",
                # Only runs when chimera.assembly is enabled, but listed
                # unconditionally like every other entry here -- the section
                # documents what this pipeline pins, not what this run used.
                "StringTie": V["stringtie"],
            },
            "TE / Gene Quantification": {
                "TEtranscripts": V["tetranscripts"],
                "TElocal": V["telocal"],
            },
            "Differential Expression & Stats": {
                "DESeq2": V["deseq2"],
                "R": V["r_base"],
            },
            # Only runs when chimera.reads is enabled, but listed
            # unconditionally like every other entry here -- see the
            # StringTie comment above.
            "Reporting": {
                "DT": V["r_dt"],
                "htmlwidgets": V["r_htmlwidgets"],
                "pandoc": V["pandoc"],
            },
        }
        with open("results/versions/rnaseq_mqc_versions.yml", "w") as fh:
            yaml.safe_dump(versions, fh, default_flow_style=False, sort_keys=False)
