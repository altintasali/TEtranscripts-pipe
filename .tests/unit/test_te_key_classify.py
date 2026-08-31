"""Feature-key parsing for the TEcount and TElocal summary sections.

These two scripts are near-twins reading different key formats, which is
exactly where they have drifted before. The keys are the only thing that
distinguishes a gene row from a TE row in a cntTable, so a parsing slip
silently reassigns whole classes of reads in the report.
"""
import pytest
import tecount_summary_mqc as tec
import telocal_summary_mqc as tel
from te_summary_common import classify


def tec_classify(key):
    return classify(key, tec.FLAVOUR.min_key_fields)


def tel_classify(key):
    return classify(key, tel.FLAVOUR.min_key_fields)


class TestTEcountKeys:
    """TEcount keys: gene_id, or gene:family:class for a TE subfamily."""

    def test_bare_gene_id_is_a_gene(self):
        assert tec_classify("ENSMUSG00000000001") == ("gene", None, None)

    def test_subfamily_key_yields_family_and_class(self):
        assert tec_classify("L1PA2:L1:LINE") == ("te", "L1", "LINE")

    def test_unrecognised_class_is_bucketed_not_dropped(self):
        kind, family, cls = tec_classify("Foo:Bar:NotARealClass")
        assert (kind, family) == ("te", "Bar")
        assert cls == "unknown"

    @pytest.mark.parametrize("cls", ["LINE", "SINE", "LTR", "DNA", "RC"])
    def test_every_declared_class_survives(self, cls):
        assert tec_classify(f"fam:fam:{cls}")[2] == cls

    def test_empty_field_is_not_a_te(self):
        # "a::b" is malformed; treating it as a TE would invent a family
        assert tec_classify("a::b") == ("gene", None, None)

    def test_two_field_key_has_no_class_of_its_own(self):
        # the last field is the family, so there is no class to read
        assert tec_classify("L1PA2:L1") == ("te", "L1PA2", "unknown")


class TestTElocalKeys:
    """TElocal keys: transcript_id:gene_id:family_id:class_id."""

    def test_bare_gene_id_is_a_gene(self):
        assert tel_classify("ENSMUSG00000000001") == ("gene", None, None)

    def test_locus_key_takes_last_two_fields(self):
        key = "chr1:564318:564741(L1PA2:+):L1PA2:L1:LINE/L1"
        kind, family, cls = tel_classify(key)
        assert kind == "te"
        assert family == "L1"
        assert cls == "unknown"  # "LINE/L1" is not one of TE_CLASSES

    def test_four_field_key(self):
        assert tel_classify("L1PA2_dup1:L1PA2:L1:LINE") == ("te", "L1", "LINE")

    def test_two_field_key_is_a_gene_for_telocal(self):
        # TElocal always has >= 3 fields for a TE; 2 means something else
        assert tel_classify("a:b") == ("gene", None, None)


def test_both_scripts_agree_on_a_shared_shape():
    """A 3-field gene:family:class key means the same thing to both."""
    assert tec_classify("x:L1:LINE") == tel_classify("x:L1:LINE")
