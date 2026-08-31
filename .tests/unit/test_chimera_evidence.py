"""The merged gene-TE evidence table.

Its contract is the pipeline's central claim about chimeras: it reports
evidence and scores nothing. The two signals deliberately excluded -- read
depth and TE-locus expression -- must never produce a flag, because both were
measured to be misleading (depth is artifact-inflated; TE expression is
anti-correlated with the splice motif).

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
    """The confidence tier was removed; it must not come back."""
    header, _ = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t2\t10\n")
    assert "confidence_tier" not in header
    assert "evidence" in header and "n_evidence" in header


def test_read_depth_earns_no_flag(tmp_path):
    """999 reads on a single sample with no motif is exactly the artifact shape."""
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tno\ttrans\tno\t1\t999\n")
    assert rows[0]["evidence"] == "."
    assert rows[0]["n_evidence"] == "0"
    assert rows[0]["junction_reads"] == "999"  # still reported


def test_te_expression_earns_no_flag(tmp_path):
    """telocal_active is anti-evidence -- reported, never counted."""
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tno\tte_initiated\tyes\t1\t5\n")
    assert rows[0]["telocal_active"] == "yes"
    assert rows[0]["evidence"] == "."


def test_canonical_and_replicates_flag(tmp_path):
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t3\t5\n")
    assert set(rows[0]["evidence"].split(",")) == {"canonical", "multi_sample"}
    assert rows[0]["n_evidence"] == "2"


def test_both_screens_flag_and_event_aggregation(tmp_path):
    _, rows = run_evidence(
        tmp_path,
        "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_terminated\tno\t2\t50\n"
        "j2\tG\tT\tsf\tfam\tLINE\tno\tte_exonized\tno\t1\t10\n",
        "MSTRG.1\tT\tsf\tfam\tLINE\tG\tyes\tte_terminated\n",
    )
    assert len(rows) == 1, "events for one pair must collapse to one row"
    row = rows[0]
    assert row["found_by"] == "both"
    assert row["junction_events"] == "2"
    assert row["junction_reads"] == "60", "reads sum across events"
    assert row["junction_canonical"] == "yes", "canonical on any event flags the pair"
    assert "both_screens" in row["evidence"]
    assert "assembly_strand_match" in row["evidence"]


def test_without_assembly_every_pair_is_junction_only(tmp_path):
    _, rows = run_evidence(tmp_path, "j1\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n")
    assert rows[0]["found_by"] == "junction"
    assert "both_screens" not in rows[0]["evidence"]


def test_sort_is_deterministic_and_unweighted(tmp_path):
    """n_evidence desc, then gene, then te -- a count, not a ranking."""
    _, rows = run_evidence(
        tmp_path,
        "j1\tZed\tT1\tsf\tfam\tLINE\tyes\tte_initiated\tno\t3\t5\n"   # 2 flags
        "j2\tAlpha\tT2\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"  # 1 flag
        "j3\tBeta\tT3\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t900\n"  # 1 flag
    )
    assert [r["gene_id"] for r in rows] == ["Zed", "Alpha", "Beta"], (
        "ties break alphabetically, not by read depth"
    )


def test_pairs_with_no_gene_or_te_are_dropped(tmp_path):
    _, rows = run_evidence(
        tmp_path,
        "j1\t.\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
        "j2\tG\t.\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
        "j3\tG\tT\tsf\tfam\tLINE\tyes\tte_initiated\tno\t1\t5\n"
    )
    assert len(rows) == 1
