#!/usr/bin/env Rscript
# Standalone, self-contained interactive HTML explorer over the FULL gene-TE
# chimera candidate catalogue (results/chimera/candidates.tsv.gz). MultiQC's
# own "Candidates" table (chimera_candidates_table_mqc.py) caps at top_n rows
# because embedding the full catalogue there runs multiqc_report.html from a
# few MB to tens of MB and slows MultiQC's own build by an order of
# magnitude (measured directly: 2.3 MB -> 83 MB, 1.4s -> 15s build, at
# real-cohort scale of ~31k rows). This is the "every row, sortable,
# searchable, filterable, no MultiQC, no server" companion.
#
# DT gives per-column search boxes for free (filter = "top"): a text box for
# character columns, a two-handle numeric range slider for numeric columns,
# and (because the low-cardinality columns below are coerced to factor()) a
# dropdown for those -- all with zero hand-written JS. htmlwidgets::
# saveWidget(selfcontained = TRUE) bundles the JS/CSS/data into ONE file via
# pandoc, so the result opens in any browser with no server and no internet
# connection needed to view it.
#
# Same non-ranking stance as the MultiQC table (chimera_candidates_table_mqc.py,
# guards 36/50): every evidence column is shown exactly as chimera_evidence.py
# wrote it, and NO combined/weighted "confidence" column is computed here --
# ever. The table opens sorted by Evidence types only because DT needs some
# initial order to open with; it is a COUNT of flags, not a score, and any
# column header re-sorts on click.
#
# Two columns are new here and don't exist in candidates.tsv.gz: "Gene locus"
# and "TE locus", a ready-to-paste IGV coordinate ("chr:start-end") for every
# row. The existing per-sample/per-transcript IGV BED tracks
# (chimera_reads_igv_bed / chimera_assembly_igv_bed) are keyed on a
# breakpoint-coordinate string or a StringTie transcript_id -- neither
# matches a candidate's gene_id/te_id, so a biologist could not search IGV by
# candidate name even with those tracks turned on (they default to off).
# results/reference/genes.bed and te.bed, by contrast, are keyed on exactly
# gene_id/te_id and are always written whenever any chimera screen is
# enabled -- joining against them gives a locus for every candidate with no
# extra config toggle.
#
# Usage:
#   Rscript chimera_candidates_explorer.R \
#       candidates.tsv.gz gene_id_to_name.tsv.gz genes.bed te.bed out.html
suppressMessages({
    library(DT)
    library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    stop(paste(
        "usage: chimera_candidates_explorer.R candidates.tsv.gz",
        "gene_id_to_name.tsv.gz genes.bed te.bed out.html"
    ))
}
candidates_path <- args[1]
gene_names_path <- args[2]
genes_bed_path <- args[3]
te_bed_path <- args[4]
out_html <- args[5]

# "." (or blank) -> NA, an empty DT cell -- not 0. chimera_evidence.py uses
# "." for telocal_count specifically to keep "TElocal never ran" (NA here)
# distinguishable from "ran, found nothing" (0). Coercing "." to 0 instead
# would silently erase that distinction.
int_or_na <- function(x) {
    x[x %in% c(".", "")] <- NA
    suppressWarnings(as.integer(x))
}

# Ready-to-paste IGV locus strings from a BED file keyed on column 4
# (genes.bed/te.bed: chrom, start, end, id, score, strand, ...; no header
# row -- see annotation_to_bed.py). BED start is 0-based half-open; IGV's
# locus box is 1-based inclusive, hence the +1 on start only.
read_bed_loci <- function(path) {
    bed <- read.delim(path, header = FALSE, colClasses = "character",
                       check.names = FALSE)
    setNames(
        sprintf("%s:%d-%s", bed[[1]], as.integer(bed[[2]]) + 1L, bed[[3]]),
        bed[[4]]
    )
}

candidates <- read.delim(candidates_path, colClasses = "character",
                          check.names = FALSE)

symbols <- tryCatch({
    s <- read.delim(gene_names_path, colClasses = "character",
                     check.names = FALSE)
    setNames(s$gene_name, s$gene_id)
}, error = function(e) character())
gene_label <- unname(symbols[candidates$gene_id])
unresolved <- is.na(gene_label) | gene_label %in% c(".", "")
gene_label[unresolved] <- candidates$gene_id[unresolved]

gene_loci <- read_bed_loci(genes_bed_path)
te_loci <- read_bed_loci(te_bed_path)
gene_locus <- unname(gene_loci[candidates$gene_id])
gene_locus[is.na(gene_locus)] <- "."
te_locus <- unname(te_loci[candidates$te_id])
te_locus[is.na(te_locus)] <- "."

# Header labels reuse chimera_candidates_table_mqc.py's wording exactly
# where that table already shows the column (Gene, TE insertion, Evidence
# types, Splice motif, Samples, Found by, Strand match, Reads, TE locus
# reads), so a reader moving between the two sees the same names.
df <- data.frame(
    "Gene" = gene_label,
    "gene_id" = candidates$gene_id,
    "TE insertion" = candidates$te_id,
    "Gene locus" = gene_locus,
    "TE locus" = te_locus,
    "TE subfamily" = candidates$te_subfamily,
    "TE family" = candidates$te_family,
    "TE class" = factor(candidates$te_class),
    "Found by" = factor(candidates$found_by),
    "Evidence flags" = candidates$evidence,
    "Evidence types" = int_or_na(candidates$n_evidence),
    "Splice motif" = factor(candidates$junction_canonical),
    "Samples" = int_or_na(candidates$junction_max_samples),
    "Junction events" = int_or_na(candidates$junction_events),
    "Reads" = int_or_na(candidates$junction_reads),
    "Junction chimera types" = candidates$junction_chimera_types,
    "TElocal active" = factor(candidates$telocal_active),
    "TE locus reads" = int_or_na(candidates$telocal_count),
    "Assembly transcripts" = int_or_na(candidates$assembly_transcripts),
    "Assembly chimera types" = candidates$assembly_chimera_types,
    "Strand match" = factor(candidates$assembly_strand_match),
    "Assembly transcript IDs" = candidates$assembly_transcript_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
# Row order is left exactly as chimera_evidence.py wrote it
# (-n_evidence, gene_id, te_id) -- not re-sorted here, so this table and
# candidates.tsv.gz always agree on order, same rationale as
# chimera_candidates_table_mqc.py's "take the head rather than re-sort".

evidence_col <- which(colnames(df) == "Evidence types") - 1L

widget <- DT::datatable(
    df,
    rownames = FALSE,
    filter = "top",
    caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left;",
        htmltools::HTML(sprintf(
            "<strong>%d</strong> gene-TE chimera candidate pairs from <code>%s</code>. ",
            nrow(df), candidates_path
        )),
        htmltools::strong("No combined score or ranking is computed here."),
        " Click a header to sort; use the boxes/sliders under the headers to filter. ",
        "Gene locus / TE locus are ready to paste into IGV's locus box."
    ),
    options = list(
        pageLength = 25,
        lengthMenu = list(c(25, 50, 100, 500, -1), c("25", "50", "100", "500", "All")),
        scrollX = TRUE,
        order = list(list(evidence_col, "desc"))
    )
)

htmlwidgets::saveWidget(widget, out_html, selfcontained = TRUE,
                         title = "Gene-TE chimera candidates")
cat(sprintf("chimera candidates explorer: %d rows -> %s\n", nrow(df), out_html))
