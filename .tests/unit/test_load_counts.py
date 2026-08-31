"""cntTable parsing in the two summary scripts.

TEcount and TElocal both emit a two-column `feature<TAB>count` table, and the
two scripts parse it with near-identical code. They had drifted: TEcount used
int(float(v)) and TElocal used int(v), so a count written as "12.0" -- which
TEtranscripts does emit, since its EM assigns fractional reads before
rounding -- crashed the TElocal path and not the TEcount one.

Both are pinned to the tolerant behaviour here.
"""
import pytest
from te_summary_common import load_counts

# Both sections now share one loader; the parametrisation is kept so the test
# still reads as "this must hold for both paths", which is what regressed.
LOADERS = [
    pytest.param(load_counts, id="tecount"),
    pytest.param(load_counts, id="telocal"),
]


def write_table(tmp_path, body, name="s.cntTable"):
    p = tmp_path / name
    p.write_text("gene/TE\tsample\n" + body)
    return p


@pytest.mark.parametrize("load", LOADERS)
def test_reads_integer_counts(load, tmp_path):
    p = write_table(tmp_path, "ENSG1\t10\nL1:L1:LINE\t5\n")
    assert load(p) == {"ENSG1": 10, "L1:L1:LINE": 5}


@pytest.mark.parametrize("load", LOADERS)
def test_float_formatted_counts_do_not_crash(load, tmp_path):
    """The regression: TEcount's EM can write "12.0" rather than "12"."""
    p = write_table(tmp_path, "ENSG1\t12.0\nENSG2\t0.0\n")
    assert load(p) == {"ENSG1": 12, "ENSG2": 0}


@pytest.mark.parametrize("load", LOADERS)
def test_header_is_skipped(load, tmp_path):
    p = write_table(tmp_path, "ENSG1\t1\n")
    assert "gene/TE" not in load(p)


@pytest.mark.parametrize("load", LOADERS)
def test_blank_lines_are_ignored(load, tmp_path):
    p = write_table(tmp_path, "ENSG1\t1\n\nENSG2\t2\n\n")
    assert load(p) == {"ENSG1": 1, "ENSG2": 2}


@pytest.mark.parametrize("load", LOADERS)
def test_empty_table_is_not_an_error(load, tmp_path):
    p = write_table(tmp_path, "")
    assert load(p) == {}
