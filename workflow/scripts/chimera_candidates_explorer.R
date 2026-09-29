#!/usr/bin/env Rscript
# Standalone, self-contained interactive HTML explorer over the FULL gene-TE
# chimera candidate catalogue (results/chimera/candidates.tsv.gz). MultiQC's
# own "Candidates" table (chimera_candidates_table_mqc.py) caps at
# --max-rows because embedding the full catalogue there runs multiqc_report.html from a
# few MB to tens of MB and slows MultiQC's own build by an order of
# magnitude (measured directly: 2.3 MB -> 83 MB, 1.4s -> 15s build, at
# real-cohort scale of ~31k rows; a real 84-sample cohort produced 454,593
# candidate rows, at which point R alone peaked at 12.3 GiB just building
# this file -- see resources.yaml). This is the "every row, sortable,
# searchable, filterable, no MultiQC, no server" companion.
#
# DT gives per-column search boxes for free (filter = "top"): a text box for
# character columns, a two-handle numeric range slider for numeric columns,
# and (because the low-cardinality columns below are coerced to factor()) a
# dropdown for those -- all with zero hand-written JS. htmlwidgets::
# saveWidget(selfcontained = TRUE) bundles the JS/CSS/data into ONE file via
# pandoc, so the result opens in any browser with no server and no internet
# connection needed to view it. Copy/CSV/Excel export buttons (DT's Buttons
# extension) and the tooltip/font tweaks below are all bundled the same way
# -- no CDN, nothing fetched at view time.
#
# Same non-ranking stance as the MultiQC table (chimera_candidates_table_mqc.py,
# guards 36/50): every evidence column is shown exactly as chimera_evidence.py
# wrote it, and NO combined/weighted "confidence" column is computed here --
# ever. The table opens sorted by Screens desc, then Gene asc, then TE
# insertion asc -- the same tie-break the MultiQC table's own defaultsort
# uses, and (unlike the Screen evidence count / Corroboration count columns
# it used to open sorted on) all three are visible by default, so what the
# table opens sorted by matches what a reader can actually see. Screens is
# a COUNT, not a score, and any column header re-sorts on click.
#
# Two columns are new here and don't exist in candidates.tsv.gz: "Gene locus"
# and "TE locus", a ready-to-paste IGV coordinate ("chr:start-end") for every
# row. The existing per-sample/per-transcript IGV BED tracks
# (chimera_chimeric_reads_igv_bed / chimera_assembly_igv_bed) are keyed on a
# breakpoint-coordinate string or a StringTie transcript_id -- neither
# matches a candidate's gene_id/te_id, so a biologist could not search IGV by
# candidate name even with those tracks turned on (they default to off).
# results/reference/genes.bed and te.bed, by contrast, are keyed on exactly
# gene_id/te_id and are always written whenever any chimera screen is
# enabled -- joining against them gives a locus for every candidate with no
# extra config toggle.
#
# Two more columns join in a real cohort-total read count, when the
# corresponding screen ran: "TElocal reads" and "Assembly reads", from
# pre-summed lookup tables built by
# chimera_candidates_matrix_totals.py (see that script and the rule in
# chimera_chimeric_reads.smk for why results/telocal/counts_matrix.tsv.gz and
# results/chimera/assembly/counts_matrix.tsv.gz are never loaded directly
# here -- both can be genome/cohort-scale).
#
# Column headers reuse chimera_candidates_table_mqc.py's wording exactly
# where that table shows the same column (kept in sync by eye -- if you
# rename a header here that also appears there, rename it there too), so a
# reader moving between the two sees the same names: CR/SJ/Assembly block
# prefixes match found_by's own cr/sj/assembly tokens, and Screen evidence
# count / Corroboration count are hidden by default here for the same
# reason that table dropped them -- every flag behind them is its own
# column now, so the counts are redundant, not deleted.
#
# Usage:
#   Rscript chimera_candidates_explorer.R \
#       candidates.tsv.gz gene_id_to_name.tsv.gz genes.bed te.bed \
#       telocal_totals.tsv assembly_totals.tsv out.html
# telocal_totals.tsv / assembly_totals.tsv may be header-only (0 data rows)
# when the corresponding screen didn't run -- every row's cohort-total join
# then resolves to NA (rendered as an empty DT cell), same convention as
# every other "not measured" cell here.
suppressMessages({
    library(DT)
    library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop(paste(
        "usage: chimera_candidates_explorer.R candidates.tsv.gz",
        "gene_id_to_name.tsv.gz genes.bed te.bed",
        "telocal_totals.tsv assembly_totals.tsv out.html"
    ))
}
candidates_path <- args[1]
gene_names_path <- args[2]
genes_bed_path <- args[3]
te_bed_path <- args[4]
telocal_totals_path <- args[5]
assembly_totals_path <- args[6]
out_html <- args[7]

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

# chimera_candidates_matrix_totals.py's output: key <TAB> total. Empty (or
# header-only, when the screen was off) reads as a zero-length named vector,
# so every lookup below naturally resolves to NA.
read_totals <- function(path) {
    d <- read.delim(path, colClasses = "character", check.names = FALSE)
    if (nrow(d) == 0) {
        return(setNames(numeric(0), character(0)))
    }
    setNames(as.numeric(d$total), d$key)
}

# assembly_transcript_ids is comma-joined (a candidate can be backed by
# several assembled transcripts) -- sum whichever of those transcript_ids
# have a total, NA (not 0) when none do (covers both "assembly didn't run"
# and "ids present but none matched", the same "not measured" convention
# used throughout this table).
sum_assembly_totals <- function(ids_str, totals) {
    vapply(ids_str, function(s) {
        if (is.na(s) || s %in% c(".", "")) {
            return(NA_real_)
        }
        vals <- totals[strsplit(s, ",", fixed = TRUE)[[1]]]
        if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
    }, numeric(1), USE.NAMES = FALSE)
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

telocal_totals <- read_totals(telocal_totals_path)
assembly_totals <- read_totals(assembly_totals_path)
telocal_cohort_total <- unname(telocal_totals[candidates$telocal_locus])
assembly_cohort_total <- sum_assembly_totals(candidates$assembly_transcript_ids,
                                              assembly_totals)

# multi_sample mirrored exactly from chimera_evidence.py's own rule
# (max(cr_max_samples, sj_max_samples) > 1) via the corroboration column it
# already wrote, not recomputed from the raw sample counts here -- same
# derivation chimera_candidates_table_mqc.py uses for its own Replicated
# column.
replicated <- factor(ifelse(
    grepl("(^|,)multi_sample(,|$)", candidates$corroboration), "yes", "no"
))

# Columns are grouped by block, in the same order and vocabulary as
# chimera_candidates_table_mqc.py: Pair, Overview, CR, SJ, Assembly,
# Support -- with this table's own explorer-only extras (loci, event
# counts, per-screen TE type, transcript IDs, flag strings) folded into
# whichever block they belong to.
df <- data.frame(
    "Gene" = gene_label,
    "gene_id" = candidates$gene_id,
    "TE insertion" = candidates$te_id,
    "Gene locus" = gene_locus,
    "TE locus" = te_locus,
    "TE subfamily" = candidates$te_subfamily,
    "TE family" = candidates$te_family,
    "TE class" = factor(candidates$te_class),
    "Screens" = int_or_na(candidates$n_screens),
    "Found by" = factor(candidates$found_by),
    "Screen evidence flags" = candidates$screen_evidence,
    "Screen evidence count" = int_or_na(candidates$n_screen_evidence),
    "Corroboration flags" = candidates$corroboration,
    "Corroboration count" = int_or_na(candidates$n_corroboration),
    "CR motif" = factor(candidates$cr_canonical),
    "CR samples" = int_or_na(candidates$cr_max_samples),
    "CR reads" = int_or_na(candidates$cr_reads),
    "CR events" = int_or_na(candidates$cr_events),
    "CR TE type" = candidates$cr_chimera_types,
    "SJ motif" = factor(candidates$sj_canonical),
    "SJ samples" = int_or_na(candidates$sj_max_samples),
    "SJ reads" = int_or_na(candidates$sj_reads),
    "SJ events" = int_or_na(candidates$sj_events),
    "SJ TE type" = candidates$sj_chimera_types,
    "Assembly strand" = factor(candidates$assembly_strand_match),
    "Assembly transcripts" = int_or_na(candidates$assembly_transcripts),
    "Assembly reads" = assembly_cohort_total,
    "Assembly TE type" = candidates$assembly_chimera_types,
    "Assembly transcript IDs" = candidates$assembly_transcript_ids,
    "Replicated" = replicated,
    "TElocal active" = factor(candidates$telocal_active),
    "TElocal reads" = telocal_cohort_total,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
# Row order is left exactly as chimera_evidence.py wrote it
# (-n_screens, -n_screen_evidence, -n_corroboration, gene_id, te_id) -- not
# re-sorted here, so this table and candidates.tsv.gz always agree on
# order. chimera_candidates_table_mqc.py's own top-N selection now sorts on
# Screens only (see that script), because it shows a capped subset and
# needs to disclose how that cap chooses rows; this table shows every row,
# so there is no such choice to make or disclose.

# One description per column, IN THE SAME ORDER as data.frame() above --
# rendered as a hover tooltip on the header (see `sketch` below). Reused
# verbatim from chimera_candidates_table_mqc.py's own COLUMNS descriptions
# for every column that table also shows, so the two views agree on wording
# there too, not just on header text.
descriptions <- c(
    "Gene symbol where the reference GTF provides one, otherwise the gene_id.",
    "Stable reference gene ID (GTF).",
    paste("The individual TE copy (transcript_id in the TE GTF), not the",
          "subfamily. Joins against TElocal rows."),
    "Genomic coordinates ready to paste into IGV (genes.bed).",
    "Genomic coordinates ready to paste into IGV (te.bed).",
    "TE annotation field from the curated TE GTF.",
    "TE annotation field from the curated TE GTF.",
    "TE annotation field from the curated TE GTF.",
    paste("How many of the 3 independent detection screens found this pair",
          "(1-3). Informational, derived from Found by."),
    paste("Which screens called it, \"+\"-joined: cr (chimeric-reads",
          "screen, STAR), assembly (StringTie), sj (SJ.out.tab screen,",
          "STAR) -- e.g. \"cr+sj\". Agreement measured near its chance",
          "rate. The CR / SJ / Assembly column-block headers below match",
          "these same three tokens."),
    paste("Screen-bound quality signals this pair carries -- see Screen",
          "evidence count. Each one can only be set if its own screen",
          "found the pair. Hidden by default: every flag here is now its",
          "own column below."),
    paste("How many screen-bound quality flags this pair carries (0-3).",
          "Can never exceed Screens (each flag needs its own screen to",
          "have found the pair). A count, not a score -- the flags are",
          "unweighted. Hidden by default, same reason as the flags column."),
    paste("Cross-cutting corroboration signals this pair carries -- see",
          "Corroboration count. Unlike Screen evidence flags, these do NOT",
          "require any particular screen. Hidden by default: every flag",
          "here is now its own column below."),
    paste("How many cross-cutting corroboration flags this pair carries",
          "(0-2: replicate support, TElocal expression). Deliberately NOT",
          "bounded by Screens -- a pair found by only one screen can",
          "still carry both. Hidden by default, same reason as the flags",
          "column."),
    paste("A recognised splice motif on at least one chimeric-junction",
          "read (STAR, chimeric-reads screen). The best artifact",
          "discriminator available."),
    "Most samples any one chimeric junction for this pair was seen in (STAR, chimeric-reads screen).",
    paste("Chimeric reads supporting this pair (STAR, chimeric-reads",
          "screen), summed across every sample that saw any of this",
          "pair's junction events -- a real cohort total, not a",
          "per-sample figure. The metric most inflated by artifacts --",
          "shown last in this block on purpose."),
    "Distinct chimeric junction events backing this pair (STAR, chimeric-reads screen).",
    paste("TE-chimera class(es) seen across this pair's junction events",
          "(STAR, chimeric-reads screen): te_initiated, te_terminated,",
          "te_exonized, or antisense_to_gene when a stranded library shows",
          "the transcript on the gene's opposite strand",
          "(see classify_chimera_chimeric_reads.py); \".\"",
          "when not classifiable (e.g. trans events on different",
          "chromosomes)."),
    paste("A recognised splice motif on at least one SJ.out.tab junction",
          "(chimera.splice_junctions screen) -- kept separate from CR",
          "motif above: a structurally independent measurement (STAR's",
          "normal-splice motif call, not the chimeric-junction one).",
          "\".\" when that screen is disabled."),
    paste("Most samples any one SJ.out.tab junction for this pair was seen",
          "in (chimera.splice_junctions screen)."),
    paste("STAR-reported unique reads across this pair's SJ.out.tab",
          "junctions (chimera.splice_junctions screen), summed across",
          "every sample that saw any of them -- a real cohort total, same",
          "caveat as CR reads above."),
    "Distinct SJ.out.tab junction events backing this pair (chimera.splice_junctions screen).",
    paste("TE-chimera class(es) seen across this pair's SJ.out.tab junction",
          "events (chimera.splice_junctions screen): te_initiated,",
          "te_terminated, te_exonized, or antisense_to_gene when the",
          "junction runs on the gene's opposite strand (see",
          "classify_chimera_splice_junctions.py); \".\" when not",
          "classifiable or that screen is disabled."),
    "The assembled transcript's strand agrees with the gene's (StringTie, assembly screen).",
    "Number of StringTie-assembled transcripts classified as this gene-TE chimera (StringTie, assembly screen).",
    paste("Sum of this pair's assembled transcript(s) estimated read count",
          "across every sample -- a real cohort total, same caveat as CR",
          "reads / SJ reads above",
          "(results/chimera/assembly/counts_matrix.tsv.gz; StringTie)."),
    paste("TE-chimera class(es) seen across this pair's assembled",
          "transcripts (StringTie, assembly screen): te_initiated,",
          "te_terminated, te_exonized, or antisense_to_gene when the",
          "transcript runs on the gene's opposite strand",
          "(see classify_chimera_assembly.py);",
          "\".\" when not classifiable (e.g. trans events on different",
          "chromosomes)."),
    "StringTie transcript_id(s) backing this pair.",
    paste("Seen in more than one sample by either the CR or SJ screen (the",
          "greater of CR samples / SJ samples is above 1). Corroboration,",
          "not screen-bound evidence -- can fire from a single screen",
          "alone. Mirrors the multi_sample flag in Corroboration flags",
          "above."),
    "Whether TElocal called this TE locus expressed in at least one sample.",
    paste("Sum of this TE locus's TElocal read count across every sample",
          "(results/telocal/counts_matrix.tsv.gz). Counted as",
          "corroboration when nonzero. Blank means TElocal did not run.")
)
stopifnot(length(descriptions) == ncol(df))

screens_col <- which(colnames(df) == "Screens") - 1L
gene_col <- which(colnames(df) == "Gene") - 1L
te_col <- which(colnames(df) == "TE insertion") - 1L
# Hidden by default (still present, searchable, exportable) -- reduces
# initial layout/render cost at high row/column counts without dropping any
# data. Computed from column NAME, not a hardcoded position, so this can't
# silently drift if a column is added/reordered above. "CR TE type" /
# "SJ TE type" / "Assembly TE type" (te_initiated/te_terminated/te_exonized)
# stay visible -- these are the TE-chimera classes readers come here looking
# for, unlike gene_id (redundant with "Gene" next to it) and the four
# Screen evidence/Corroboration flag/count columns (redundant with the CR
# motif / SJ motif / Assembly strand / Replicated / TElocal reads columns
# now shown individually -- same reason chimera_candidates_table_mqc.py
# dropped them from its own table entirely rather than hiding them).
hidden_cols <- which(colnames(df) %in% c(
    "gene_id", "Screen evidence flags", "Screen evidence count",
    "Corroboration flags", "Corroboration count"
)) - 1L

header_titles <- descriptions
sketch <- htmltools::withTags(table(
    class = "display",
    thead(
        tr(lapply(seq_along(colnames(df)), function(i) {
            th(colnames(df)[i], title = header_titles[i])
        }))
    )
))

widget <- DT::datatable(
    df,
    container = sketch,
    rownames = FALSE,
    filter = "top",
    extensions = "Buttons",
    caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left;",
        htmltools::HTML(sprintf(
            "<strong>%d</strong> gene-TE chimera candidate pairs from <code>%s</code>. ",
            nrow(df), candidates_path
        )),
        htmltools::strong("No combined score or ranking is computed here."),
        " Click a header to sort (hover a header for what it means); use the",
        " boxes/sliders under the headers to filter. Gene locus / TE locus",
        " are ready to paste into IGV's locus box.",
        htmltools::tags$br(), htmltools::tags$br(),
        htmltools::em(
            "Acronyms: TE = transposable element. In Found by / Evidence",
            "flags, cr = chimeric-reads screen (STAR), sj = SJ.out.tab",
            "screen (STAR), assembly = StringTie screen."
        )
    ),
    options = list(
        pageLength = 25,
        lengthMenu = list(c(25, 50, 100, 500, -1), c("25", "50", "100", "500", "All")),
        scrollX = TRUE,
        order = list(list(screens_col, "desc"), list(gene_col, "asc"),
                     list(te_col, "asc")),
        # Large-table performance: deferRender skips per-row work until a
        # row is actually displayed; autoWidth off skips DataTables'
        # automatic column-width measurement pass across every row/column.
        # Both are additive, zero behavior change to sort/search/filter, no
        # new dependency (r-dt already ships both). Neither shrinks the
        # embedded data payload -- see the script-level comment above on
        # why that's close to a hard floor for this single-file design.
        deferRender = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(visible = FALSE, targets = hidden_cols)),
        dom = "Blfrtip",
        buttons = list("copy", "csv", "excel")
    )
)

# Modern sans-serif stack (Helvetica/Arial, with the OS-native UI fonts
# preferred where available) applied to the whole page, not just the table
# -- injected as a plain <style> tag rather than a DT `class`/theme option,
# so it survives saveWidget(selfcontained = TRUE) as ordinary embedded CSS.
widget <- htmlwidgets::appendContent(widget, htmltools::tags$style(htmltools::HTML(
    paste0(
        "body, table.dataTable, .dataTables_wrapper { font-family: ",
        "-apple-system, 'Segoe UI', 'Helvetica Neue', Helvetica, Arial, ",
        "sans-serif; }"
    )
)))

htmlwidgets::saveWidget(widget, out_html, selfcontained = TRUE,
                         title = "Gene-TE chimera candidates")
cat(sprintf("chimera candidates explorer: %d rows -> %s\n", nrow(df), out_html))
