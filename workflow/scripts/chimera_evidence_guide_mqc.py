#!/usr/bin/env python3
"""The report's guide to reading the gene-TE chimera evidence -- and nothing
more than a guide.

This pipeline does not rank chimera candidates. It used to: a four-tier
confidence ladder, then three competing per-screen top-N tables, both
removed for having no validated weighting (commits d927c8f, 0d04e43). The
Candidates table above (chimera_candidates_table_mqc.py) is the sortable,
unranked replacement -- every pair, every evidence column, ordered by
Screens only, re-orderable by clicking any header. This section explains
what each of that table's columns is worth and how to read it; it never
lists a candidate pair itself, so it cannot become a second, differently
ordered ranking.

Project measurements that used to be hard-coded here -- numbers from one
specific run, printed into every OTHER run's report regardless of species
or cohort, some of them already stale by the time a reader saw them -- now
live in docs/chimera-evidence.md instead, with the cohort, the restriction
each measurement was made under, and what it does and doesn't show. This
guide states only what is true for any run: how each signal is defined and
what its known failure modes are.

Emits two MultiQC custom-content documents:

  chimera_evidence_guide_mqc.json        html   how to weigh each signal
  chimera_evidence_composition_mqc.json  bar    pairs carrying each signal
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

PARENT_ID = "chimera"
PARENT_NAME = "Chimera"

# Order mirrors the Candidates table's own column order (CR motif, SJ
# motif, Assembly strand, Replicated, TElocal reads) -- keep this in sync
# with chimera_candidates_table_mqc.py's COLUMNS if that table's columns
# ever change. Drives the composition bar's category order (sort_samples:
# False below, since MultiQC otherwise alphabetises bar categories) and
# this file's own signal-table row order (see signals()).
FLAGS = [
    ("cr_canonical", "CR motif"),
    ("sj_canonical", "SJ motif"),
    ("assembly_strand_match", "Assembly strand"),
    ("multi_sample", "Replicated"),
    ("telocal_expressed", "TElocal expressed"),
]

WEIGHT_STYLE = {
    "strong": ("#1a7f5a", "Best discriminator"),
    "mixed": ("#8a6d15", "Read with care"),
    "not_validated": ("#8a3a54", "Not validated"),
    "not-evidence": ("#777", "Not evidence"),
}


def _sj_motif_row(sj_require_canonical):
    """The SJ motif signal row. Its "how to read it" text depends on this
    run's chimera.splice_junctions.require_canonical: when true, the flag
    is nearly guaranteed by construction and says little on its own; when
    false, it is a genuine per-junction measurement. Rendered conditionally
    so the guide states the correct case for THIS run rather than hedging
    both at once."""
    if sj_require_canonical:
        how = (
            "chimera.splice_junctions.require_canonical is true in this "
            "run, so nearly every SJ-screen pair already carries this "
            "flag by construction &mdash; it discriminates little within "
            "this screen's own output, and is mainly useful for "
            "cross-checking against the other two screens (see "
            "Screens / Found by above)."
        )
    else:
        how = (
            "chimera.splice_junctions.require_canonical is false in this "
            "run, so this flag is not guaranteed by construction and "
            "should discriminate normally within this screen's own "
            "output -- check its distribution in your own "
            "candidates.tsv.gz."
        )
    return (
        "SJ motif",
        "STAR (SJ.out.tab)",
        "A recognised splice motif on at least one normal splice junction "
        "from the splice_junctions screen (<code>sj_canonical</code>) "
        "&mdash; kept separate from CR motif above: a structurally "
        "independent measurement (chimeric-junction typing vs. SJ.out.tab "
        "motif).",
        how,
        "mixed",
    )


def _replicated_row(star_two_pass):
    """The Replicated signal row. Under star.two_pass: cohort, junctions
    from every sample's first pass are inserted into every sample's index,
    so a junction is easier to find again in the other samples: SJ samples
    and Replicated are then not fully independent detections. Rendered for
    THIS run's setting, like the SJ motif row; a generic caveat when the
    setting is not passed."""
    how = (
        "A sequence-driven template switch or ligation artifact recurs "
        "across libraries too, so recurrence alone does not separate a "
        "real chimera from a reproducible one."
    )
    if star_two_pass == "cohort":
        how += (
            " This run used <code>star.two_pass: cohort</code>: junctions "
            "found in any sample's first pass are added to every sample's "
            "index, which makes a junction easier to detect again in the "
            "other samples. SJ samples and Replicated are therefore not "
            "fully independent detections here &mdash; the effect is small "
            "but only ever upward."
        )
    elif star_two_pass in ("per_sample", "none"):
        how += (
            f" This run used <code>star.two_pass: {star_two_pass}</code>, "
            "so each sample's junctions were detected independently of the "
            "other samples."
        )
    else:
        how += (
            " With <code>star.two_pass: cohort</code>, junctions pooled "
            "from every sample's first pass are added to every sample's "
            "index, so SJ samples and Replicated are not fully independent "
            "detections."
        )
    return (
        "Replicated",
        "STAR (chimeric junctions or SJ.out.tab)",
        "Seen in more than one sample, from the greater of the "
        "CR samples / SJ samples columns (<code>multi_sample</code>) "
        "&mdash; corroboration, not screen-bound evidence, so it can "
        "fire from a single screen alone.",
        how,
        "mixed",
    )


def signals(sj_require_canonical, star_two_pass=None):
    """(label, source tool, what the signal is, how to read it, standing)
    for the guide table, in the SAME order as FLAGS / the Candidates
    table's own columns (TE orientation first: its TE-vs-gene block sits
    before Screens there) -- one deliberate exception: Read depth spans
    both the CR and SJ blocks and is not evidence at all, so it goes last
    rather than being split into two rows."""
    return [
        (
            "TE orientation",
            "Annotation (genes.bed, te.bed)",
            "Whether the TE copy lies on the same strand as the gene "
            "(<code>te_orientation</code>: sense / antisense), shown with "
            "where it sits relative to the gene (<code>te_position</code>: "
            "upstream, intronic, exonic, downstream) and how far away "
            "(<code>te_gene_distance_bp</code>). From the annotation, so it "
            "needs no stranded library.",
            "Read it only together with TE position: each position has its "
            "own background mix of orientations, so a sense or antisense "
            "TE means little without knowing where it sits. An LTR "
            "promoter driving a gene should be sense to it. Intronic TEs "
            "have their own baseline &mdash; sense-oriented L1s, for "
            "example, are depleted from introns genome-wide &mdash; and "
            "exonized SINEs are often antisense. Compare groups of pairs "
            "against that background; a single pair's orientation proves "
            "nothing.",
            "mixed",
        ),
        (
            "Screens / Found by",
            "STAR + StringTie + STAR (SJ.out.tab)",
            "How many of the 3 independent detection screens found this "
            "pair (<code>n_screens</code> / Screens, 1-3), and which "
            "ones (<code>found_by</code> / Found by).",
            "The three screens have different blind spots, so agreement "
            "between them should in principle be strong support. How "
            "much it actually adds has not been established. When the "
            "SJ screen requires canonical junctions "
            "(<code>chimera.splice_junctions.require_canonical</code>), "
            "agreement with it partly selects for the motif by "
            "construction rather than confirming it independently.",
            "not_validated",
        ),
        (
            "Status (novel / annotated / antisense)",
            "All three screens",
            "Which kinds of chimera support the pair "
            "(<code>chimera_status</code> / Status), \"+\"-joined: "
            "<b>novel</b> (a new TE-initiated, TE-terminated or TE-exonized "
            "transcript), <b>annotated</b> (a TE-driven transcript the "
            "reference annotation already has: a TE promoter, a TE in a "
            "terminal exon, or an annotated splice into a TE-derived exon) "
            "and <b>antisense</b> (a transcript joining the TE to the "
            "gene's exon on the opposite strand).",
            "Annotated chimeras are real but known, so they are not new "
            "findings. Antisense chimeras cannot make the gene's mRNA but "
            "may regulate the gene. The SJ and assembly screens take the "
            "strand from the splice motif, so they call antisense on any "
            "library; the chimeric-reads screen uses the read strand and "
            "can only call antisense with a stranded library, so on "
            "unstranded data some of its novel calls may be antisense. A "
            "pair with several kinds sums their reads in its per-screen "
            "columns.",
            "mixed",
        ),
        (
            "CR motif",
            "STAR (chimeric junctions)",
            "A recognised splice motif on at least one chimeric-junction "
            "read (<code>cr_canonical</code>).",
            "Real splice junctions are almost always canonical; "
            "template-switching, ligation and PCR chimeras usually carry "
            "no motif at all. Compare rates between groups of pairs, not "
            "against an absolute number -- a low overall rate on its own "
            "is normal.",
            "strong",
        ),
        (
            "CR max anchor",
            "STAR (chimeric junctions)",
            "The best chimeric read's shorter segment, in aligned bases "
            "(<code>cr_max_anchor</code>), from the chimeric-reads calls "
            "only.",
            "A breakpoint that no read anchors well on both sides is easy "
            "to produce by mis-mapping or a chance alignment of a short "
            "fragment. A short maximum anchor makes that possible for this "
            "pair. A signal to read alongside the others, never a filter "
            "or a score.",
            "mixed",
        ),
        _sj_motif_row(sj_require_canonical),
        (
            "SJ unique fraction / overhang",
            "STAR (SJ.out.tab)",
            "How much of a pair's SJ-screen support maps uniquely "
            "(<code>sj_unique_fraction</code>: unique reads / all reads) and "
            "the longest anchor any read had across the junction "
            "(<code>sj_max_overhang</code>), from the SJ chimera calls only.",
            "Reads from young TE families map equally well to many copies, "
            "and junctions can arise from that mis-mapping rather than real "
            "splicing. A low unique fraction or a short maximum overhang "
            "means that is possible for this pair. Both are signals to read "
            "alongside the others, never a filter or a score.",
            "mixed",
        ),
        (
            "Assembly strand",
            "StringTie (assembly)",
            "The assembled transcript's strand agrees with the gene's "
            "(<code>assembly_strand_match</code>).",
            "A consistency check on the assembly call rather than "
            "independent support. A transcript on the gene's opposite "
            "strand is antisense transcription through that gene, not a "
            "chimera of it: all three screens now type such calls "
            "<code>antisense_to_gene</code> (the chimeric-reads screen "
            "only with a stranded library) instead of "
            "<code>te_initiated</code>/<code>te_terminated</code>/"
            "<code>te_exonized</code>.",
            "mixed",
        ),
        (
            "Assembly last-exon TE distance",
            "StringTie (assembly)",
            "For a TE in an assembled transcript's last exon, how far it "
            "sits from that exon's splice acceptor "
            "(<code>assembly_te_acceptor_distance_bp</code>; 0 = the TE "
            "takes the splice).",
            "The assembly screen calls a TE-terminated transcript when "
            "any TE overlaps a new last exon, and a last exon can be a "
            "long 3' UTR. At 0 the TE supplies the splice acceptor: a "
            "TE-derived terminal exon. A short distance past the acceptor, "
            "with the TE in the transcript's own orientation, is where a "
            "TE can supply the polyadenylation signal that ends the "
            "transcript. Far from the acceptor, the TE is UTR content and "
            "the call says little about termination. Sort or filter on it; "
            "where the line falls has been measured on one cohort only "
            "(see the project's evidence notes).",
            "mixed",
        ),
        _replicated_row(star_two_pass),
        (
            "TElocal reads",
            "TElocal",
            "TElocal's read count for the TE copy itself, and whether it "
            "is called expressed in at least one sample "
            "(<code>telocal_expressed</code>) &mdash; TElocal is a "
            "fourth data source, not one of the three detection screens.",
            "An expressed TE copy makes a chimera possible but does not "
            "show one; highly expressed copies also generate more "
            "chimeric-looking artifacts simply by generating more reads.",
            "not_validated",
        ),
        (
            "Read depth (CR reads, SJ reads)",
            "STAR",
            "<strong>Not an evidence flag.</strong> Reported as "
            "<code>cr_reads</code> / CR reads and <code>sj_reads</code> "
            "/ SJ reads.",
            "The metric most inflated by artifacts &mdash; a hot PCR "
            "chimera is often the deepest event in a run. Depth never "
            "promotes a pair on its own; a high-depth row with no other "
            "flags means exactly that.",
            "not-evidence",
        ),
    ]


def load(path):
    """Rows of a TSV as dicts, keyed by header name."""
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def guide_html(sj_require_canonical, star_two_pass=None):
    rows = []
    for label, source, what, how, weight in signals(sj_require_canonical,
                                                    star_two_pass):
        colour, badge = WEIGHT_STYLE[weight]
        rows.append(
            "<tr>"
            f'<td style="white-space:nowrap;"><strong>{label}</strong></td>'
            f'<td style="white-space:nowrap;color:#555;">{source}</td>'
            f'<td style="white-space:nowrap;color:{colour};">{badge}</td>'
            f"<td>{what}</td>"
            f"<td>{how}</td>"
            "</tr>"
        )

    return f"""
<p>What each column of the <strong>Candidates</strong> table above is worth,
and how to read it &mdash; this is the reference for choosing which signal
your question needs, and how far to trust it.</p>

<table class="table" style="width:100%; font-size: 90%;">
<thead><tr>
<th>Signal</th><th>Source</th><th>Standing</th><th>What it is</th>
<th>How to read it</th>
</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>

<p style="font-size: 85%; color: #888;">Full catalogue, one row per pair
with every evidence column: <code>results/chimera/candidates.tsv.gz</code>.
Project measurements so far are recorded in
<code>docs/chimera-evidence.md</code>, and none of them has established a
weighting.</p>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True,
                    help="chimera candidates.tsv.gz (chimera_evidence.py)")
    ap.add_argument("--sj-require-canonical", required=True,
                    choices=["true", "false"],
                    help="config chimera.splice_junctions.require_canonical "
                    "-- picks the SJ motif row's conditional wording")
    ap.add_argument("--star-two-pass", default=None,
                    choices=["cohort", "per_sample", "none"],
                    help="config star.two_pass -- picks the Replicated row's "
                    "conditional wording (generic caveat when omitted)")
    ap.add_argument("--out-guide", required=True)
    ap.add_argument("--out-composition", required=True)
    args = ap.parse_args()
    sj_require_canonical = args.sj_require_canonical == "true"

    rows = list(load(args.evidence))
    n_pairs = len(rows)

    composition = {flag: 0 for flag, _ in FLAGS}
    n_no_flags = 0
    for r in rows:
        # screen_evidence and corroboration are the two independent counts
        # chimera_evidence.py splits flags into (see its module docstring):
        # the first is bounded by n_screens, the second deliberately is not.
        # This section explains all five flags together regardless of which
        # count they belong to, so both are read here.
        present = [
            f for col in ("screen_evidence", "corroboration")
            for f in r.get(col, ".").split(",") if f != "."
        ]
        if not present:
            n_no_flags += 1
        for flag in present:
            if flag in composition:
                composition[flag] += 1

    guide_doc = {
        # NOT "chimera_evidence_guide": a section id must never equal a
        # parent_id. report_section_order's first pass matches MODULE anchors,
        # and a custom-content group's anchor is its parent_id -- so a section
        # sharing that name gets picked up by the module pass, whose order
        # semantics are inverted, and the whole group renders LAST. Measured.
        #
        # Every view in this group must still use the same parent_id, or
        # MultiQC silently splits the group into two report sections.
        "id": "chimera_signal_guide",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "section_name": "How to weigh this evidence",
        "description": (
            "What each line of chimera evidence is worth, and how to "
            "read it -- click Help for the full signal-by-signal table. "
            "No ranking is produced."
        ),
        "helptext": guide_html(sj_require_canonical, args.star_two_pass),
        "plot_type": "html",
        "data": (
            "<p>Sort the <strong>Candidates</strong> table above on the "
            "signal your question needs; click <strong>Help</strong> "
            "(top right of this section) for what each column is worth "
            "and how to read it.</p>"
        ),
    }

    # Composition, not a ranking: how much of each signal the cohort produced.
    # Shown because "148 of 2,431 carry a splice motif" changes how the whole
    # table should be read, and is invisible from any individual row.
    # A bargraph whose every value is zero does not render as an empty plot --
    # MultiQC raises ValueError("No datasets to plot"), exits non-zero, and
    # NO REPORT IS WRITTEN AT ALL. A run that finds no gene-TE pairs is a
    # normal outcome (a clean library, or a species with a sparse TE
    # annotation), so it must not cost the user their entire report. Fall back
    # to an HTML section, the same way chimera_assembly_summary_mqc.py does
    # for its own empty case.
    counts_by_label = {
        **{label: composition[flag] for flag, label in FLAGS},
        "No evidence flag": n_no_flags,
    }
    if any(counts_by_label.values()):
        # One bar PER EVIDENCE TYPE, not one stacked bar split by type.
        # MultiQC stacks by default (stacking="relative"), which drew these as
        # segments of a single bar -- and that reads as a partition: mutually
        # exclusive slices summing to the cohort. It is the opposite of true.
        # A pair can carry all five flags at once, so the counts OVERLAP and
        # do not sum to anything meaningful. Giving each its own bar removes
        # the implied exclusivity.
        composition_body = {
            "plot_type": "bargraph",
            "pconfig": {
                "id": "chimera_evidence_composition_plot",
                "title": "Gene-TE pairs carrying each line of evidence",
                "ylab": "Gene-TE pairs",
                "cpswitch": False,
                "stacking": "group",
                # whole pairs: no decimals. tt_decimals is the key MultiQC
                # honours here; "format" is silently dropped.
                "tt_decimals": 0,
                # MultiQC alphabetises bar categories by default
                # (bargraph.py's sort_samples: True), which would scramble
                # FLAGS' own order above. Measured.
                "sort_samples": False,
            },
            "categories": ["Gene-TE pairs"],
            "data": {
                label: {"Gene-TE pairs": count}
                for label, count in counts_by_label.items()
            },
        }
    else:
        composition_body = {
            "plot_type": "html",
            "data": (
                "<p>No gene-TE pairs were found in this run, so there is "
                "nothing to plot here. This is not an error: it means neither "
                "screen called a gene-TE chimera. The screens are independent "
                "of each other and of TE quantification, so TEcount and "
                "TElocal results are unaffected.</p>"
            ),
        }

    composition_doc = {
        "id": "chimera_evidence_composition",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "section_name": "Evidence composition",
        "description": (
            "How many gene-TE pairs carry each line of evidence. "
            "<strong>These bars overlap and do not sum to the cohort</strong> "
            "-- a pair can carry every flag at once, so they are independent "
            "counts, not slices of a whole."
        ),
        "helptext": (
            "Sources, in the same order as the bars: CR motif from STAR "
            "chimeric junctions, SJ motif from STAR's SJ.out.tab (when "
            "chimera.splice_junctions is enabled), Assembly strand from "
            "StringTie, Replicated from STAR chimeric junctions or "
            "SJ.out.tab, TElocal expressed from TElocal. How many "
            "screens agreed on a pair is shown separately (the Screens "
            "column in the Candidates table above), not as a bar here "
            "-- it is derived from found_by, not an independent evidence "
            "flag."
        ),
        **composition_body,
    }

    for path, doc in ((args.out_guide, guide_doc),
                      (args.out_composition, composition_doc)):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open_write(path) as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")

    summary = ", ".join(f"{flag}: {composition[flag]}" for flag, _ in FLAGS)
    print(f"chimera evidence guide: {n_pairs} gene-TE pairs ({summary}, "
          f"no flag: {n_no_flags}) -> {args.out_guide}, {args.out_composition}")


if __name__ == "__main__":
    main()
