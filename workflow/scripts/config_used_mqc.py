#!/usr/bin/env python3
"""Write the resolved run configuration as a MultiQC custom-content JSON table.

Snakemake ``script:`` directive calls this with a ``snakemake`` object whose
attributes mirror the rule's config, output, and log declarations.
"""
import json
import os
import sys
import traceback

try:
    snakemake
except NameError:
    snakemake = None


def _read_version(snakemake_obj):
    """Return the pipeline version string, or 'unknown' if VERSION is missing."""
    try:
        snakefile = str(snakemake_obj.snakefile)
    except AttributeError:
        snakefile = ""
    if snakefile:
        version_path = os.path.join(
            os.path.dirname(os.path.dirname(snakefile)), "VERSION"
        )
        try:
            with open(version_path) as fh:
                return fh.read().strip()
        except OSError:
            pass
    return "unknown"


def main(smk):
    config = smk.config
    params = smk.params
    log_path = str(smk.log)

    # Ensure the log directory exists so the fallback logger below works.
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    star_extra = config.get("star", {}).get("extra", "") or "(none)"
    trim_extra = config.get("trimming", {}).get("extra", "") or "(none)"
    te_extra = config.get("tetranscripts", {}).get("extra", "") or "(none)"

    samples = params.get("_samples", [])
    sample_count = params.get("_sample_count", len(samples))
    sample_names = ", ".join(samples) if samples else ""

    sjdb_overhang = params.get("_sjdb_overhang", "auto")
    star_index = params.get("_star_index", "")
    trim_enabled = params.get("_trim_enabled", "")
    tecount_qc_enabled = params.get("_tecount_qc_enabled", "")
    tecount_qc = params.get("_tecount_qc", {})
    chimera_enabled = params.get("_chimera_enabled", "")
    telocal_enabled = params.get("_telocal_enabled", "")
    telocal_locind_auto = params.get("_telocal_locind_auto", "")
    telocal_qc_enabled = params.get("_telocal_qc_enabled", "")
    telocal_qc = params.get("_telocal_qc", {})
    keep_merged = params.get("_keep_merged_fastq", "")
    keep_trimmed = params.get("_keep_trimmed_fastq", "")
    keep_star_index = params.get("_keep_star_index", "")
    keep_telocal_index = params.get("_keep_telocal_index", "")

    rows = {
        "pipeline_version": _read_version(smk),
        "samples": f"{sample_count} ({sample_names})",
        "ref.fasta": str(config.get("ref", {}).get("fasta", "(not provided)")),
        "ref.gtf": str(config["ref"]["gtf"]),
        "ref.te_gtf": str(config["ref"]["te_gtf"]),
        "ref.sjdb_overhang": str(sjdb_overhang),
        "star.index": star_index,
        "star.extra": star_extra,
        "trimming.enabled": str(trim_enabled),
        "trimming.trim_nextseq": str(config.get("trimming", {}).get("trim_nextseq", 0)),
        "trimming.extra": trim_extra,
        "strandedness.min_fraction": str(config.get("strandedness", {}).get("min_fraction", 0.6)),
        "tetranscripts.mode": config["tetranscripts"]["mode"],
        "tetranscripts.padj": str(config["tetranscripts"]["padj"]),
        "tetranscripts.foldchange": str(config["tetranscripts"]["foldchange"]),
        "tetranscripts.minread": str(config["tetranscripts"]["minread"]),
        "tetranscripts.extra": te_extra,
        "tetranscripts.qc.enabled": str(tecount_qc_enabled),
        "tetranscripts.qc.feature_class": tecount_qc.get("feature_class", ""),
        "tetranscripts.qc.pca_transform": tecount_qc.get("pca_transform", ""),
        "chimera.reads.enabled": str(chimera_enabled),
        "chimera.reads.breakpoint_tolerance": str(config.get("chimera", {}).get("reads", {}).get("breakpoint_tolerance", 0)),
        "chimera.reads.require_canonical_junction": str(config.get("chimera", {}).get("reads", {}).get("require_canonical_junction", False)),
        "chimera.reads.qc.pca_transform": config.get("chimera", {}).get("reads", {}).get("qc", {}).get("pca_transform", "vst"),
        "chimera.assembly.enabled": str(config.get("chimera", {}).get("assembly", {}).get("enabled", False)),
        "chimera.assembly.breakpoint_tolerance": str(config.get("chimera", {}).get("assembly", {}).get("breakpoint_tolerance", 0)),
        "telocal.enabled": str(telocal_enabled),
        "telocal.locind": (
            "(auto-build from TE GTF)" if telocal_locind_auto in (True, "True", "")
            else str(config.get("telocal", {}).get("locind", ""))
        ),
        "telocal.qc.enabled": str(telocal_qc_enabled),
        "telocal.qc.feature_class": telocal_qc.get("feature_class", ""),
        "telocal.qc.pca_transform": telocal_qc.get("pca_transform", ""),
        "outputs.keep_merged_fastq": str(keep_merged),
        "outputs.keep_trimmed_fastq": str(keep_trimmed),
        "outputs.keep_star_index": str(keep_star_index),
        "outputs.keep_telocal_index": str(keep_telocal_index),
    }

    doc = {
        "id": "config_used",
        "section_name": "Configuration Used",
        "description": (
            "The resolved run configuration (config values from the "
            "config file passed with --configfile)."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "config_used_table",
            "title": "Configuration used",
            "col1_header": "Setting",
            "sort_rows": False,
        },
        "headers": {"value": {"title": "Value"}},
        "data": {k: {"value": v} for k, v in rows.items()},
    }

    # Ensure the output directory exists.
    os.makedirs(os.path.dirname(str(smk.output)), exist_ok=True)

    with open(str(smk.output), "w") as fh:
        json.dump(doc, fh, indent=2)


if __name__ == "__main__" or snakemake is not None:
    target = snakemake
    if target is None:
        print("This script must be run by Snakemake's script: directive.",
              file=sys.stderr)
        sys.exit(1)
    try:
        main(target)
    except Exception:
        # Write a full traceback to the log file so it is never empty,
        # even when Snakemake's own error capture fails.
        log_path = str(target.log)
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "a") as log_fh:
            traceback.print_exc(file=log_fh)
        raise
