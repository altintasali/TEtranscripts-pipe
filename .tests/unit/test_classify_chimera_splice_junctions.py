"""chimera_type classification for the SJ.out.tab screen.

Bug fixed here (2026): chimera_type used to be decided from where the TE's
own span sits relative to the GENE's overall genomic span (TE entirely
upstream of the gene -> te_initiated; entirely downstream -> te_terminated;
anywhere else, including squarely inside an intron -> te_exonized). That
missed the main case this screen exists to catch: a TE sitting in an
intron, acting as an alternative promoter that splices directly into a
downstream exon (donor in TE, acceptor in a gene exon -- direction
te_to_gene), used to be scored te_exonized purely because its coordinates
fall inside the gene's span, even though the junction itself says "this
transcript starts here". In mouse oocytes/early embryos, intronic
MT2/MERVL promoters are a large share of the real TE-initiated transcripts.

Fixed classification (see classify_chimera_splice_junctions.py's module
docstring for the full reasoning): decided from the junction's DIRECTION
plus the gene's own annotated EXON STRUCTURE (exons.bed), not TE-vs-gene-
span containment. chimera_type's vocabulary is unchanged
(te_initiated/te_terminated/te_exonized, matching classify_chimera_assembly
.py) so existing consumers of the column are not broken; the extra
upstream-vs-internal distinction assembly does not make is reported in a
new, separate te_initiated_detail column instead of being folded into
chimera_type.

These tests run the real CLI end to end (small BED/SJ.out.tab fixtures),
the same way test_chimera_evidence.py drives chimera_evidence.py, because
the classification logic depends on gene_exon_positions built from the
loaded exons.bed inside main() -- there is no smaller unit to call directly
without duplicating that setup.
"""
import subprocess
import sys
from pathlib import Path

SCRIPT = (Path(__file__).resolve().parents[2] / "workflow" / "scripts"
          / "classify_chimera_splice_junctions.py")


def classify(tmp_path, genes, exons, te, sj):
    (tmp_path / "genes.bed").write_text(genes)
    (tmp_path / "exons.bed").write_text(exons)
    (tmp_path / "te.bed").write_text(te)
    (tmp_path / "SJ.out.tab").write_text(sj)
    out = tmp_path / "out.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPT),
         "--sj", str(tmp_path / "SJ.out.tab"),
         "--genes", str(tmp_path / "genes.bed"),
         "--exons", str(tmp_path / "exons.bed"),
         "--te", str(tmp_path / "te.bed"),
         "--sample", "S1", "--min-unique-reads", "1",
         "--out", str(out)],
        check=True, capture_output=True,
    )
    with open(out) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return [dict(zip(header, line)) for line in lines]


# --- "+" strand ---------------------------------------------------------

def test_te_upstream_of_the_genes_own_first_exon_is_te_initiated_upstream(tmp_path):
    """Classic case: TE sits before the gene's only/first exon and splices
    straight into it -- donor in TE, acceptor at the gene's own most-5'
    annotated exon, nothing else annotated further upstream."""
    rows = classify(
        tmp_path,
        genes="chr1\t99\t6000\tXyz\t.\t+\n",
        exons="chr1\t1999\t2300\tXyz\t.\t+\n",
        te="chr1\t100\t500\tTE1\t.\t+\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t501\t1999\t1\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "te_to_gene"
    assert r["chimera_type"] == "te_initiated"
    assert r["te_initiated_detail"] == "upstream"


def test_intronic_te_promoter_skipping_an_earlier_exon_is_te_initiated_internal(tmp_path):
    """The bug's exact repro shape: TE sits inside an intron and splices
    into a downstream (non-first) exon, skipping the gene's actual first
    exon. Old logic: te_exonized (TE coordinates fall inside the gene
    span). Fixed: te_initiated, detail=internal."""
    rows = classify(
        tmp_path,
        genes="chr1\t1999\t6000\tAbc\t.\t+\n",
        exons="chr1\t1999\t2300\tAbc\t.\t+\nchr1\t4999\t5300\tAbc\t.\t+\n",
        te="chr1\t2999\t3500\tMT2_Mm_dup2\t.\t+\tERVL\tLTR\tMT2_Mm\n",
        sj="chr1\t3401\t4999\t1\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "te_to_gene"
    assert r["chimera_type"] == "te_initiated"
    assert r["te_initiated_detail"] == "internal"


def test_te_exonized_when_the_gene_continues_past_the_te(tmp_path):
    """gene_to_te (donor in a gene exon, acceptor in TE): exonized, not
    terminated, when the annotation has a further exon downstream of the
    TE -- transcription is not known to stop here."""
    rows = classify(
        tmp_path,
        genes="chr1\t99\t8000\tXyz\t.\t+\n",
        exons="chr1\t199\t500\tXyz\t.\t+\nchr1\t6999\t7300\tXyz\t.\t+\n",
        te="chr1\t2999\t3500\tTE1\t.\t+\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t501\t2999\t1\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "gene_to_te"
    assert r["chimera_type"] == "te_exonized"
    assert r["te_initiated_detail"] == "."


def test_te_terminated_when_nothing_annotated_follows_the_te(tmp_path):
    """gene_to_te with no further annotated exon past the TE -- nothing
    known follows, so this plausibly ends the transcript."""
    rows = classify(
        tmp_path,
        genes="chr1\t99\t4000\tXyz\t.\t+\n",
        exons="chr1\t199\t500\tXyz\t.\t+\n",
        te="chr1\t2999\t3500\tTE1\t.\t+\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t501\t2999\t1\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "gene_to_te"
    assert r["chimera_type"] == "te_terminated"
    assert r["te_initiated_detail"] == "."


# --- "-" strand (transcript order runs high -> low genomic coordinate) --

def test_te_upstream_minus_strand(tmp_path):
    rows = classify(
        tmp_path,
        genes="chr1\t99\t8000\tXyz\t.\t-\n",
        exons="chr1\t199\t500\tXyz\t.\t-\n",
        te="chr1\t2999\t3500\tTE1\t.\t-\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t501\t2999\t2\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "te_to_gene"
    assert r["chimera_type"] == "te_initiated"
    assert r["te_initiated_detail"] == "upstream"


def test_intronic_te_promoter_minus_strand(tmp_path):
    """Same skipped-first-exon shape as the "+" repro, mirrored: the gene's
    real first exon (highest coordinate on "-") lies beyond the TE, and the
    junction splices into a lower-coordinate, non-first exon."""
    rows = classify(
        tmp_path,
        genes="chr1\t99\t8000\tXyz\t.\t-\n",
        exons="chr1\t6999\t7300\tXyz\t.\t-\nchr1\t199\t500\tXyz\t.\t-\n",
        te="chr1\t2999\t3500\tTE1\t.\t-\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t501\t2999\t2\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "te_to_gene"
    assert r["chimera_type"] == "te_initiated"
    assert r["te_initiated_detail"] == "internal"


def test_te_exonized_minus_strand(tmp_path):
    rows = classify(
        tmp_path,
        genes="chr1\t99\t8000\tXyz\t.\t-\n",
        exons="chr1\t6999\t7300\tXyz\t.\t-\nchr1\t199\t500\tXyz\t.\t-\n",
        te="chr1\t2999\t3500\tTE1\t.\t-\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t3001\t6999\t2\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "gene_to_te"
    assert r["chimera_type"] == "te_exonized"
    assert r["te_initiated_detail"] == "."


def test_te_terminated_minus_strand(tmp_path):
    rows = classify(
        tmp_path,
        genes="chr1\t99\t8000\tXyz\t.\t-\n",
        exons="chr1\t6999\t7300\tXyz\t.\t-\n",
        te="chr1\t2999\t3500\tTE1\t.\t-\tLTR\tLTR\tTE1sub\n",
        sj="chr1\t3001\t6999\t2\t1\t0\t9\t1\t33\n",
    )
    assert len(rows) == 1
    r = rows[0]
    assert r["direction"] == "gene_to_te"
    assert r["chimera_type"] == "te_terminated"
    assert r["te_initiated_detail"] == "."


def test_chimera_type_vocabulary_matches_the_assembly_screen(tmp_path):
    """No new values leak into chimera_type itself -- cross-screen "TE
    type" comparisons in the report only make sense if both screens draw
    from the same label set. The extra upstream/internal distinction lives
    in te_initiated_detail instead."""
    rows = classify(
        tmp_path,
        genes="chr1\t1999\t6000\tAbc\t.\t+\n",
        exons="chr1\t1999\t2300\tAbc\t.\t+\nchr1\t4999\t5300\tAbc\t.\t+\n",
        te="chr1\t2999\t3500\tMT2_Mm_dup2\t.\t+\tERVL\tLTR\tMT2_Mm\n",
        sj="chr1\t3401\t4999\t1\t1\t0\t9\t1\t33\n",
    )
    assert rows[0]["chimera_type"] in ("te_initiated", "te_terminated", "te_exonized")
