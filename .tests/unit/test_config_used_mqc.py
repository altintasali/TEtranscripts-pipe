"""config_used_mqc.py's table must not silently drop a whole chimera
sub-screen (or an important toggle inside one) the way it dropped
chimera.splice_junctions entirely, and chimera.assembly.require_tss_in_te,
before this fix -- see the module's chimera.* rows.

Scope: this test walks config.schema.yaml's chimera.* subtree only, not the
whole schema. The report's config_used_mqc.py deliberately curates *every*
section (tetranscripts.qc/telocal.qc show enabled + feature_class +
pca_transform, never the qc.min_samples_present/min_total_counts/min_events
filter floors or per-rule star.*/outputs.* internals) -- that curation is an
established, repeated design choice across every screen, not a bug, and
forcing the whole schema into this table would just re-introduce the bloat
Task 3 is trying to remove. What WAS a bug: one sub-screen's block (
splice_junctions) missing outright, and one obviously report-relevant toggle
(assembly.require_tss_in_te, which changes te_initiated's definition)
missing from a block that otherwise has its siblings. ALLOWED_OMISSIONS below
is the explicit list of leaves this test does not require, grouped by why.
"""
import sys
import types
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "config_used_mqc.py"
SCHEMA = REPO_ROOT / "workflow" / "schemas" / "config.schema.yaml"

sys.path.insert(0, str(SCRIPT.parent))

# Leaves deliberately never shown in config_used_mqc.py, grouped by why --
# the same category of omission already applies uniformly to
# tetranscripts.qc/telocal.qc, not just chimera, so these are a pre-existing
# pattern, not new slack introduced by this test.
ALLOWED_OMISSIONS = {
    # STAR chimeric-alignment tuning internals: fixed, rarely-changed
    # constants, not something a report reader decides the run's validity
    # from the way they would enabled/breakpoint_tolerance/require_canonical.
    "chimera.chimeric_reads.star.segment_min",
    "chimera.chimeric_reads.star.overhang_min",
    "chimera.chimeric_reads.star.score_drop_max",
    "chimera.chimeric_reads.star.extra",
    # QC-view filter floors: same omission pattern as tetranscripts.qc /
    # telocal.qc, which also show only pca_transform, never these.
    "chimera.chimeric_reads.qc.min_samples_present",
    "chimera.chimeric_reads.qc.min_total_counts",
    "chimera.chimeric_reads.qc.min_events",
    "chimera.splice_junctions.qc.min_samples_present",
    "chimera.splice_junctions.qc.min_total_counts",
    "chimera.splice_junctions.qc.min_events",
    # Output-writing toggles (write this file or don't) -- not a modeling
    # choice that changes what a candidate means, unlike require_canonical/
    # require_tss_in_te.
    "chimera.chimeric_reads.outputs.write_igv_bed",
    "chimera.chimeric_reads.outputs.write_counts_matrix",
    "chimera.chimeric_reads.outputs.write_candidates_explorer",
    "chimera.assembly.outputs.write_igv_bed",
    "chimera.assembly.outputs.write_gene_te_chimera_counts",
    "chimera.splice_junctions.outputs.write_counts_matrix",
    # min_transcript_tpm: explicitly documented in the schema as "not
    # applied at detection time... no rule currently filters on it
    # automatically" -- a README-example threshold, not a run setting.
    "chimera.assembly.min_transcript_tpm",
}


def _schema_leaves(node, prefix=""):
    """Yield dotted paths for every scalar (non-object) property under node."""
    props = node.get("properties", {})
    for key, sub in props.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(sub, dict) and sub.get("type") == "object":
            yield from _schema_leaves(sub, path)
        else:
            yield path


def _config_used_rows():
    import config_used_mqc

    config = {
        "chimera": {
            "chimeric_reads": {"breakpoint_tolerance": 0, "require_canonical": False,
                                "qc": {"pca_transform": "vst"}},
            "assembly": {"enabled": True, "breakpoint_tolerance": 0,
                         "require_tss_in_te": True},
            "splice_junctions": {"enabled": True, "breakpoint_tolerance": 0,
                                  "min_unique_reads": 1, "require_canonical": True,
                                  "qc": {"pca_transform": "vst"}},
        },
        "ref": {"gtf": "g.gtf", "te_gtf": "te.gtf"},
        "strandedness": {"min_fraction": 0.8},
        "tetranscripts": {"mode": "multi"},
    }
    params = {
        "_samples": ["S1"], "_sample_count": 1, "_sjdb_overhang": "auto",
        "_star_index": "idx", "_trim_enabled": True,
        "_tecount_qc_enabled": True, "_tecount_qc": {},
        "_chimera_chimeric_reads_enabled": True,
        "_telocal_enabled": True, "_telocal_locind_auto": True,
        "_telocal_qc_enabled": True, "_telocal_qc": {},
        "_keep_merged_fastq": True, "_keep_trimmed_fastq": True,
        "_keep_star_index": True, "_keep_telocal_index": True,
    }
    tmp_log = "/tmp/test_config_used_mqc.log"
    tmp_out = "/tmp/test_config_used_mqc.json"
    smk = types.SimpleNamespace(config=config, params=params, log=tmp_log, output=tmp_out)
    config_used_mqc.main(smk)
    import json
    with open(tmp_out) as fh:
        doc = json.load(fh)
    return set(doc["data"].keys())


def test_every_chimera_schema_leaf_is_shown_or_explicitly_omitted():
    import config_used_mqc

    with open(SCHEMA) as fh:
        schema = yaml.safe_load(fh)
    chimera_node = schema["properties"]["chimera"]
    leaves = {f"chimera.{p}" for p in _schema_leaves(chimera_node)}

    rows = _config_used_rows()
    shown = {config_used_mqc.row_label(leaf) for leaf in leaves} & rows
    shown_dotted = {leaf for leaf in leaves if config_used_mqc.row_label(leaf) in rows}
    missing = leaves - shown_dotted - ALLOWED_OMISSIONS
    assert not missing, (
        f"chimera.* schema leaves missing from config_used_mqc.py's table "
        f"and not in ALLOWED_OMISSIONS: {sorted(missing)}"
    )
    assert shown  # sanity: the transform above actually matched something

    # Guard the allowlist itself: every entry must still be a real schema
    # leaf (catches a stale entry left behind after a schema rename).
    stale = ALLOWED_OMISSIONS - leaves
    assert not stale, f"ALLOWED_OMISSIONS entries no longer in the schema: {sorted(stale)}"


def test_assembly_and_splice_junctions_report_their_defining_toggle():
    """Pin the two fields this fix specifically added, by name."""
    import config_used_mqc

    rows = _config_used_rows()
    assert config_used_mqc.row_label("chimera.assembly.require_tss_in_te") in rows
    assert config_used_mqc.row_label("chimera.splice_junctions.enabled") in rows
    assert config_used_mqc.row_label("chimera.splice_junctions.require_canonical") in rows


def test_row_labels_survive_multiqcs_sample_name_cleaning():
    """The actual bug: MultiQC's custom-content loader runs every row key
    through its generic sample-name cleaner (base_module.clean_s_name),
    unconditionally, using config.fn_clean_exts/fn_clean_trim -- MultiQC's
    own ~150+ built-in patterns PLUS this repo's multiqc_config.yaml
    extra_fn_clean_exts appended on top. Reproduces that exact merge and
    checks every row_label() this script emits survives intact. Skips
    cleanly if multiqc isn't importable in this environment (not a hard
    dependency of the unit test suite)."""
    multiqc = pytest.importorskip("multiqc")
    from multiqc import config as mqc_config
    from multiqc.base_module import BaseMultiqcModule

    with open(REPO_ROOT / "workflow" / "default-config" / "multiqc_config.yaml") as fh:
        mqc_yaml = yaml.safe_load(fh)
    extra = mqc_yaml.get("extra_fn_clean_exts", [])
    mqc_config.fn_clean_exts = mqc_config.fn_clean_exts + extra

    class _Fake(BaseMultiqcModule):
        def __init__(self):
            pass

    bm = _Fake()
    rows = _config_used_rows()
    mangled = []
    for label in rows:
        f = {"root": "", "fn": label, "sp_key": None}
        cleaned = bm.clean_s_name(label, f)
        if cleaned != label:
            mangled.append((label, cleaned))
    assert not mangled, f"row labels mangled by MultiQC's sample-name cleaning: {mangled}"
