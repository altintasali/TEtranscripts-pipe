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


def row_label(key):
    """Turn a dotted config-path key into this table's row label/id.

    See the "data" field's comment in main() for why this exists: MultiQC's
    custom-content loader runs every row key through its generic sample-name
    cleaner regardless of whether the key is actually a sample name, and
    that cleaner's ~150+ built-in patterns (plus this repo's own
    extra_fn_clean_exts) are almost all anchored on a literal "." or "_" --
    both of which are all over a dotted config path. Neutralizing them
    avoids the whole collision surface instead of chasing individual
    patterns. Exposed as its own function so
    .tests/unit/test_config_used_mqc.py can verify against the exact
    transform actually used, not a re-implementation of it.
    """
    return key.replace(".", " > ").replace("_", " ")


def _read_version():
    """Return the pipeline version string, or 'unknown' if VERSION is missing.

    Resolved from this script's own __file__ (workflow/scripts/ -> workflow
    -> repo root), not from the `snakemake` script object: Snakemake's
    `script:` directive does not expose a `.snakefile` attribute (only
    input/output/params/wildcards/threads/resources/log/config/rule/
    scriptdir/bench_iteration are), so the previous attempt to read
    `snakemake_obj.snakefile` always raised AttributeError and silently
    returned "unknown" -- verified: the field showed "unknown" in every
    generated config_used_mqc.json. workflow/scripts/tetranscripts-pipe and
    .tests/check_version_sync.py already solve this exact problem correctly
    via __file__; this follows the same pattern (and, via __file__, still
    resolves through a `workflow/` symlink to a shared checkout -- see the
    tetranscripts-pipe CLI's own docstring for that deployment shape).
    """
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        with open(os.path.join(repo_root, "VERSION")) as fh:
            return fh.read().strip()
    except OSError:
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
    chimera_enabled = params.get("_chimera_chimeric_reads_enabled", "")
    telocal_enabled = params.get("_telocal_enabled", "")
    telocal_locind_auto = params.get("_telocal_locind_auto", "")
    telocal_qc_enabled = params.get("_telocal_qc_enabled", "")
    telocal_qc = params.get("_telocal_qc", {})
    keep_merged = params.get("_keep_merged_fastq", "")
    keep_trimmed = params.get("_keep_trimmed_fastq", "")
    keep_star_index = params.get("_keep_star_index", "")
    keep_telocal_index = params.get("_keep_telocal_index", "")

    rows = {
        "pipeline_version": _read_version(),
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
        "strandedness.min_fraction": str(config["strandedness"]["min_fraction"]),
        "strandedness.balanced_max": str(config.get("strandedness", {}).get("balanced_max", 0.55)),
        "tetranscripts.mode": config["tetranscripts"]["mode"],
        "tetranscripts.extra": te_extra,
        "tetranscripts.qc.enabled": str(tecount_qc_enabled),
        "tetranscripts.qc.feature_class": tecount_qc.get("feature_class", ""),
        "tetranscripts.qc.pca_transform": tecount_qc.get("pca_transform", ""),
        "chimera.chimeric_reads.enabled": str(chimera_enabled),
        "chimera.chimeric_reads.breakpoint_tolerance": str(config.get("chimera", {}).get("chimeric_reads", {}).get("breakpoint_tolerance", 0)),
        "chimera.chimeric_reads.require_canonical": str(config.get("chimera", {}).get("chimeric_reads", {}).get("require_canonical", False)),
        "chimera.chimeric_reads.qc.enabled": str(config.get("chimera", {}).get("chimeric_reads", {}).get("qc", {}).get("enabled", False)),
        "chimera.chimeric_reads.qc.pca_transform": config.get("chimera", {}).get("chimeric_reads", {}).get("qc", {}).get("pca_transform", "vst"),
        "chimera.assembly.enabled": str(config.get("chimera", {}).get("assembly", {}).get("enabled", False)),
        "chimera.assembly.breakpoint_tolerance": str(config.get("chimera", {}).get("assembly", {}).get("breakpoint_tolerance", 0)),
        "chimera.assembly.require_tss_in_te": str(config.get("chimera", {}).get("assembly", {}).get("require_tss_in_te", True)),
        "chimera.splice_junctions.enabled": str(config.get("chimera", {}).get("splice_junctions", {}).get("enabled", False)),
        "chimera.splice_junctions.breakpoint_tolerance": str(config.get("chimera", {}).get("splice_junctions", {}).get("breakpoint_tolerance", 0)),
        "chimera.splice_junctions.min_unique_reads": str(config.get("chimera", {}).get("splice_junctions", {}).get("min_unique_reads", 1)),
        "chimera.splice_junctions.require_canonical": str(config.get("chimera", {}).get("splice_junctions", {}).get("require_canonical", True)),
        "chimera.splice_junctions.qc.enabled": str(config.get("chimera", {}).get("splice_junctions", {}).get("qc", {}).get("enabled", False)),
        "chimera.splice_junctions.qc.pca_transform": config.get("chimera", {}).get("splice_junctions", {}).get("qc", {}).get("pca_transform", "vst"),
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
        # BUG FIXED 2026: row_label() neutralizes "." and "_" so MultiQC's
        # generic sample-name cleaner (applied unconditionally to every
        # dict-shaped custom-content row key, sample name or not -- see
        # row_label()'s docstring) can't truncate these labels. Verified:
        # not only this repo's own extra_fn_clean_exts caused this --
        # "trimming.trim_nextseq" collided with MultiQC's OWN default
        # ".trim" pattern (for Trim Galore log filenames) and was truncated
        # to "trimming" even before this repo's settings applied.
        "data": {row_label(k): {"value": v} for k, v in rows.items()},
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
