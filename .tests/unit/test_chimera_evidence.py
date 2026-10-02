"""The merged gene-TE evidence table.

Its contract is the pipeline's central claim about chimeras: it reports
evidence and scores nothing. Read depth is the one signal deliberately
excluded from ever producing a flag, because it was measured to be
misleading (depth is artifact-inflated). TE-locus expression (telocal_active)
IS counted as a flag despite a mixed real-cohort measurement -- see
chimera_evidence.py's module docstring for the full reasoning.

The flag logic lives inside main(), so this drives the CLI the way the rule
does and reads the table back.
"""
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "chimera_evidence.py"

JUNCTION_HEADER = ("event_id\tgene_id\tte_id\tte_subfamily\tte_family\tte_class\t"
                   "canonical\tchimera_type\ttelocal_active\tn_samples\ttotal_reads\n")
ASSEMBLY_HEADER = ("transcript_id\tte_id\tte_subfamily\tte_family\tte_class\t"
                   "matched_gene_id\tstrand_match\tchimera_type\n")


def gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def run_evidence(tmp_path, junction_rows, assembly_rows=None):
    j = gz(tmp_path / "j.tsv.gz", JUNCTION_HEADER + junction_rows)
    out = tmp_path / "out.tsv.gz"
    cmd = [sys.executable, str(SCRIPT), "--junction", str(j), "--out", str(out)]
    if assembly_rows is not None:
        a = gz(tmp_path / "a.tsv.gz", ASSEMBLY_HEADER + assembly_rows)
        cmd += ["--assembly", str(a)]
    subprocess.run(cmd, check=True, capture_output=True)
    with gzip.open(out, "rt") as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh]
    header, rows = lines[0], lines[1:]
    return header, [dict(zip(header, r)) for r in rows]


def test_no_scoring_column_exists(tmp_path):
    """The confidence tier was removed; it must not come back. evidence/
    n_evidence were themselves replaced by the screen_evidence/corroboration
    split (see chimera_evidence.py's module docstring) -- they must not come
    back either, since a single combined count mixed screen-bound flags
    (bounded by n_screens) with cross-cutting ones (not bounded by it)."""
    header, _ = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t2\t10\n")
    assert "confidence_tier" not in header
    assert "evidence" not in header and "n_evidence" not in header
    assert "screen_evidence" in header and "n_screen_evidence" in header
    assert "corroboration" in header and "n_corroboration" in header


def test_read_depth_earns_no_flag(tmp_path):
    """999 reads on a single sample with no motif is exactly the artifact shape."""
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tno\tte_initiated\tno\t1\t999\n")
    assert rows[0]["screen_evidence"] == "."
    assert rows[0]["n_screen_evidence"] == "0"
    assert rows[0]["corroboration"] == "."
    assert rows[0]["n_corroboration"] == "0"
    assert rows[0]["cr_reads"] == "999"  # still reported


def test_te_expression_earns_a_corroboration_flag(tmp_path):
    """telocal_active is an unresolved-but-counted signal -- see
    chimera_evidence.py's module docstring for why it stays a flag despite
    the mixed real-cohort measurement. It is corroboration, not screen
    evidence: TElocal is a fourth data source, not one of the three
    detection screens."""
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tno\tte_initiated\tyes\t1\t5\n")
    assert rows[0]["telocal_active"] == "yes"
    assert rows[0]["screen_evidence"] == "."
    assert rows[0]["corroboration"] == "telocal_expressed"
    assert rows[0]["n_corroboration"] == "1"


def test_canonical_and_replicates_split_across_the_two_counts(tmp_path):
    """cr_canonical is screen-bound (needs the cr screen to have found the
    pair); multi_sample is cross-cutting (fires from one screen alone) --
    they must land in different columns, not one combined count."""
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t3\t5\n")
    assert rows[0]["screen_evidence"] == "cr_canonical"
    assert rows[0]["n_screen_evidence"] == "1"
    assert rows[0]["corroboration"] == "multi_sample"
    assert rows[0]["n_corroboration"] == "1"


def test_n_screens_and_event_aggregation(tmp_path):
    """n_screens replaces the old both_screens/all_three_screens flags: it is
    a plain count derived from found_by, not summed into either evidence
    count below (see chimera_evidence.py's module docstring)."""
    _, rows = run_evidence(
        tmp_path,
        "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_terminated\tno\t2\t50\n"
        "j2\tG\tT\tsf\tfam\tLINE\tno\tte_exonized\tno\t1\t10\n",
        "MSTRG.1\tT\tsf\tfam\tLINE\tG\tyes\tte_terminated\n",
    )
    assert len(rows) == 1, "events for one pair must collapse to one row"
    row = rows[0]
    assert row["found_by"] == "cr+assembly"
    assert row["n_screens"] == "2", "reads + assembly found this pair"
    assert row["cr_events"] == "2"
    assert row["cr_reads"] == "60", "reads sum across events"
    assert row["cr_canonical"] == "yes", "canonical on any event flags the pair"
    assert set(row["screen_evidence"].split(",")) == {"cr_canonical", "assembly_strand_match"}
    assert row["n_screen_evidence"] == "2", "bounded by n_screens: both screens found it"
    assert "both_screens" not in row["screen_evidence"], "removed: duplicated n_screens/found_by"
    assert "all_three_screens" not in row["screen_evidence"]


def test_without_assembly_every_pair_is_junction_only(tmp_path):
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n")
    assert rows[0]["found_by"] == "cr"
    assert rows[0]["n_screens"] == "1"
    assert "both_screens" not in rows[0]["screen_evidence"]


def test_sort_is_deterministic_and_unweighted(tmp_path):
    """(n_screens, n_screen_evidence, n_corroboration) desc, then gene, then
    te -- three counts, not a ranking. All three rows here are cr-only
    (n_screens=1, tied) and all carry a canonical motif (n_screen_evidence=1,
    tied), so the tie-break that actually separates them is
    n_corroboration (Zed's n_samples=3 earns multi_sample; Alpha/Beta's
    n_samples=1 does not), then alphabetical."""
    _, rows = run_evidence(
        tmp_path,
        "j1\tZed\tT1\tsf\tfam\tLINE\tyes\tte_initiated\tno\t3\t5\n"    # n_corroboration=1
        "j2\tAlpha\tT2\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"  # n_corroboration=0
        "j3\tBeta\tT3\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t900\n"  # n_corroboration=0
    )
    assert [r["gene_id"] for r in rows] == ["Zed", "Alpha", "Beta"], (
        "ties break alphabetically, not by read depth"
    )
    assert [r["n_corroboration"] for r in rows] == ["1", "0", "0"]


def test_pairs_with_no_gene_or_te_are_dropped(tmp_path):
    _, rows = run_evidence(
        tmp_path,
        "j1\t.\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
        "j2\tG\t.\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
        "j3\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
    )
    assert len(rows) == 1
