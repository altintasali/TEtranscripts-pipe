"""The interval search behind every gene/TE annotation lookup.

classify_chimera_junctions.overlapping() decides which genes and TEs a chimeric
breakpoint lands in, so an off-by-one or a missed long feature here silently
mislabels calls rather than failing. Two real bugs in this area are in the
history ("breakpoints were resolved one base off", "te.bed collapsed every TE
subfamily into one genome-spanning interval"), which is why the awkward cases
-- nested, adjacent, and long-preceding features -- are pinned explicitly.

Coordinates are half-open [start, end), matching BED.
"""
import classify_chimera_junctions as pcj
import pytest


def track(*features):
    """Build the {chrom: (feats, max_end)} structure load_bed produces."""
    feats = sorted((s, e, extras) for s, e, extras in features)
    max_end, running = [], 0
    for _, e, _ in feats:
        running = max(running, e)
        max_end.append(running)
    return {"chr1": (feats, max_end)}


def names(hits):
    return sorted(h[2][0] for h in hits)


class TestBasicOverlap:
    def test_hit_inside_a_feature(self):
        t = track((100, 200, ("A",)))
        assert names(pcj.overlapping(t, "chr1", 150, 151)) == ["A"]

    def test_unknown_chromosome_is_empty_not_an_error(self):
        t = track((100, 200, ("A",)))
        assert pcj.overlapping(t, "chrZ", 150, 151) == []

    def test_query_entirely_before_and_after(self):
        t = track((100, 200, ("A",)))
        assert pcj.overlapping(t, "chr1", 0, 50) == []
        assert pcj.overlapping(t, "chr1", 300, 400) == []


class TestHalfOpenBoundaries:
    """[start, end) -- the exact place an off-by-one hides."""

    def test_query_touching_feature_start_overlaps(self):
        t = track((100, 200, ("A",)))
        assert names(pcj.overlapping(t, "chr1", 99, 101)) == ["A"]

    def test_query_ending_exactly_at_feature_start_does_not_overlap(self):
        t = track((100, 200, ("A",)))
        assert pcj.overlapping(t, "chr1", 50, 100) == []

    def test_query_starting_exactly_at_feature_end_does_not_overlap(self):
        t = track((100, 200, ("A",)))
        assert pcj.overlapping(t, "chr1", 200, 250) == []

    def test_last_base_of_a_feature_overlaps(self):
        t = track((100, 200, ("A",)))
        assert names(pcj.overlapping(t, "chr1", 199, 200)) == ["A"]


class TestAwkwardLayouts:
    def test_nested_feature_is_found(self):
        t = track((100, 900, ("outer",)), (400, 500, ("inner",)))
        assert names(pcj.overlapping(t, "chr1", 450, 451)) == ["inner", "outer"]

    def test_long_feature_starting_well_before_the_query(self):
        """The case a naive bisect on `start` misses."""
        t = track((0, 10_000, ("long",)), (100, 200, ("short",)))
        assert names(pcj.overlapping(t, "chr1", 5000, 5001)) == ["long"]

    def test_short_feature_sandwiched_between_long_ones(self):
        t = track((0, 10_000, ("long",)), (5000, 5100, ("mid",)),
                  (9000, 20_000, ("long2",)))
        assert names(pcj.overlapping(t, "chr1", 5050, 5051)) == ["long", "mid"]

    def test_adjacent_features_do_not_bleed(self):
        t = track((100, 200, ("A",)), (200, 300, ("B",)))
        assert names(pcj.overlapping(t, "chr1", 150, 160)) == ["A"]
        assert names(pcj.overlapping(t, "chr1", 250, 260)) == ["B"]

    def test_query_spanning_several_features(self):
        t = track((100, 200, ("A",)), (200, 300, ("B",)), (300, 400, ("C",)))
        assert names(pcj.overlapping(t, "chr1", 150, 350)) == ["A", "B", "C"]


class TestStrandHelper:
    @pytest.mark.parametrize("given,expected",
                             [("+", "-"), ("-", "+"), (".", "."), ("", "."),
                              ("?", ".")])
    def test_opp(self, given, expected):
        assert pcj.opp(given) == expected
