"""Fisher's exact test and BH correction behind the splice-motif enrichment.

scipy is not in this workflow's environment and is far too heavy a dependency
for one test, so chimera_reads_qc_mqc implements both. That makes them OUR
correctness problem: a wrong p-value here would not fail anything, it would
just publish a wrong number in the report next to a real one.

Every expected value below was produced by R -- `fisher.test(matrix(...))` and
`p.adjust(..., method = "BH")` -- and is pinned to 10 significant figures.
"""
import chimera_reads_qc_mqc as jq
import pytest

# (a, b, c, d) -> R fisher.test(matrix(c(a, c, b, d), nrow = 2))$p.value
R_FISHER = [
    ((3, 1, 1, 3), 0.4857142857),          # tea tasting, the textbook case
    ((10, 3, 2, 15), 0.0005367241191),
    ((1, 9, 11, 3), 0.002759456185),
    ((100, 200, 150, 400), 0.07011451857),
    ((5, 0, 0, 5), 0.007936507937),        # zero cells
    ((1647, 29250, 900, 20000), 1.119619402e-07),  # cohort-sized counts
    ((12, 5, 8, 20), 0.01224926264),
]


@pytest.mark.parametrize("table,expected", R_FISHER)
def test_matches_r_fisher_test(table, expected):
    assert jq.fisher_exact_two_sided(*table) == pytest.approx(expected, rel=1e-9)


def test_p_is_a_probability():
    for table, _ in R_FISHER:
        assert 0.0 <= jq.fisher_exact_two_sided(*table) <= 1.0


def test_symmetric_under_transpose():
    """Swapping the two rows cannot change a two-sided p."""
    for a, b, c, d in [(10, 3, 2, 15), (100, 200, 150, 400), (12, 5, 8, 20)]:
        assert jq.fisher_exact_two_sided(a, b, c, d) == pytest.approx(
            jq.fisher_exact_two_sided(c, d, a, b), rel=1e-12)


def test_empty_margin_is_not_significant():
    """A margin with no observations cannot support any claim."""
    assert jq.fisher_exact_two_sided(0, 0, 5, 5) == 1.0
    assert jq.fisher_exact_two_sided(0, 5, 0, 5) == 1.0


def test_bh_matches_r_p_adjust():
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    expected = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1,
                0.1057142857, 0.216, 0.216, 0.216]
    assert jq.benjamini_hochberg(p) == pytest.approx(expected, rel=1e-9)


def test_bh_is_monotone_and_order_preserving():
    """q must follow p's ordering and never fall below it."""
    p = [0.3, 0.01, 0.2, 0.04, 0.9]
    q = jq.benjamini_hochberg(p)
    assert all(qi >= pi for qi, pi in zip(q, p))
    assert [q[i] for i in sorted(range(len(p)), key=p.__getitem__)] == sorted(q)


def test_bh_handles_empty():
    assert jq.benjamini_hochberg([]) == []
