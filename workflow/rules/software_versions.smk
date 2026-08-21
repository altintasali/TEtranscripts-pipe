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

        # Every pinned tool, under the pipeline's group name -- MultiQC
        # renders the top-level YAML key as the group header in its Software
        # Versions section. Tools whose versions MultiQC auto-detects from
        # logs (STAR, FastQC, samtools...) are listed here too, so the
        # section stays complete even when a log omits them.
        versions = {
            "STAR": V["star"],
            "samtools": V["samtools"],
            "Trim Galore!": V["trim_galore"],
            "FastQC": V["fastqc"],
            "RSeQC": V["rseqc"],
            "MultiQC": V["multiqc"],
            "TEtranscripts": V["tetranscripts"],
            "TElocal": V["telocal"],
            "DESeq2": V["deseq2"],
            "pheatmap": V["pheatmap"],
            "R": V["r_base"],
            "UCSC tools": f"{V['ucsc_gtftogenepred']}",
        }
        with open("results/versions/rnaseq_mqc_versions.yml", "w") as fh:
            yaml.safe_dump(
                {"TEtranscripts-pipe": versions}, fh, default_flow_style=False
            )
