"""GTF attribute parsing, duplicated in two scripts that must agree.

Both annotation_to_bed.py and gene_name_lookup.py pull gene_id/gene_name out
of a GTF attribute column. A parsing slip here is how "te.bed collapsed every
TE subfamily into one genome-spanning interval" happened: the wrong attribute
was used as the interval key.
"""
import annotation_to_bed as a2b
import gene_name_lookup as gnl

GENCODE = 'gene_id "ENSMUSG01"; gene_name "Actb"; gene_type "protein_coding";'
REPEATMASKER = 'gene_id "L1PA2"; transcript_id "L1PA2_dup1"; family_id "L1"; class_id "LINE";'


class TestAnnotationToBed:
    def test_extracts_quoted_values(self):
        attrs = a2b.parse_attrs(GENCODE)
        assert attrs["gene_id"] == "ENSMUSG01"
        assert attrs["gene_name"] == "Actb"

    def test_repeatmasker_style_keeps_transcript_id_distinct(self):
        """transcript_id is the individual insertion; gene_id is the subfamily.
        Conflating them is what collapsed te.bed."""
        attrs = a2b.parse_attrs(REPEATMASKER)
        assert attrs["gene_id"] == "L1PA2"
        assert attrs["transcript_id"] == "L1PA2_dup1"
        assert attrs["gene_id"] != attrs["transcript_id"]

    def test_missing_attribute_is_absent_not_empty(self):
        assert "gene_name" not in a2b.parse_attrs('gene_id "X";')

    def test_empty_attribute_column(self):
        assert a2b.parse_attrs("") == {}


class TestGeneNameLookup:
    def test_extracts_quoted_values(self):
        attrs = gnl.parse_attrs(GENCODE)
        assert attrs["gene_id"] == "ENSMUSG01"
        assert attrs["gene_name"] == "Actb"

    def test_value_containing_spaces_survives(self):
        attrs = gnl.parse_attrs('gene_id "X"; gene_name "some long name";')
        assert attrs["gene_name"] == "some long name"

    def test_agrees_with_annotation_to_bed_on_shared_keys(self):
        a, g = a2b.parse_attrs(GENCODE), gnl.parse_attrs(GENCODE)
        for key in ("gene_id", "gene_name"):
            assert a[key] == g[key]
