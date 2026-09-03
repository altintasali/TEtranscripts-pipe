#!/usr/bin/env Rscript
# Sample-QC: normalize a counts matrix (chimera junctions, TEcount features
# or TElocal loci) and produce the PCA + sample-distance views shipped by
# chimera_reads_qc.smk / tecount_qc.smk / telocal.smk.
#
# The first two arguments select the mode and the view being served:
#   view   "chimera", "tecount" or "telocal" -- namespaces the MultiQC
#          custom-content ids/titles so all views can render in one report
#          without colliding.
# Two modes, selected by the script's argument vector:
#   --transform view counts.tsv samples.csv transform min_samples_present \
#                 min_total_counts out_matrix.tsv
#       Apply the chosen transformation (vst / rlog / log2) to the counts
#       matrix (feature x sample, as written by chimera_reads_counts.py /
#       tecount_counts.py) and write the transformed matrix for the plot
#       rule to read. vst/rlog use DESeq2's blind normalization
#       (independent of sample labels, so the QC view can't be overfit);
#       log2 is log2(x + 1) without DESeq2. Features seen in fewer than
#       min_samples_present samples, or with fewer than min_total_counts
#       supporting reads overall, are dropped for this view only.
#   --plots view transformed.tsv samples.csv min_events transform \
#                 out_pca_mqc.json out_heatmap_mqc.json
#       PCA scatter + sample-to-sample Euclidean-distance heatmap of the
#       transformed matrix, written as MultiQC custom-content JSON documents
#       (the `_mqc.json` suffix makes multiqc_report.html render them
#       interactively). PCA points are colored by the sample sheet's
#       "condition" column; samples without a condition form one group.
#
# QC filters (chimera.qc / tetranscripts.qc in the config) apply ONLY to this
# view: the full catalog/counts matrix is never reduced.
# Libraries are loaded per-mode (DESeq2 for the transform); the plots mode
# needs no extra packages -- the JSON payloads are small and built by hand.
suppressMessages(library(DESeq2))

# Per-view naming: the MultiQC custom-content ids, section names, plot titles
# and descriptions are namespaced per view so the chimera, TEcount and TElocal
# QC views can all render inside one multiqc_report.html without colliding.
# Section names carry their own grouping now. The chimera views share ONE
# parent group with every other chimera section (parent = "chimera"), because
# MultiQC has no third heading level -- so a bare "PCA" there would sit among a
# dozen sibling sections with nothing saying which screen it belongs to. The
# tecount/telocal views keep bare names: they have a parent group to themselves.
#
# `id` namespaces the emitted doc ids and must stay unique per view; `parent`
# is the MultiQC group and must match the Python emitters' parent_id exactly.
VIEWS <- list(
    # NB: the list KEY stays "chimera" -- it is the CLI selector passed by
    # chimera_reads_qc.smk (--transform chimera / --plots chimera).
    chimera = list(
        id = "chimera_reads",
        parent = "chimera",
        label = "Chimera",
        section_prefix = "Reads - ",
        noun_plural = "chimeric events",
        noun_singular = "event"
    ),
    assembly = list(
        id = "chimera_assembly",
        parent = "chimera",
        label = "Chimera",
        section_prefix = "Assembly - ",
        noun_plural = "assembled chimeric transcripts",
        noun_singular = "transcript"
    ),
    tecount = list(
        id = "tecount",
        parent = "tecount",
        label = "TEcount",
        section_prefix = "",
        noun_plural = "features",
        noun_singular = "feature"
    ),
    telocal = list(
        id = "telocal",
        parent = "telocal",
        label = "TElocal",
        section_prefix = "",
        noun_plural = "loci",
        noun_singular = "locus"
    )
)

view_params <- function(view) {
    if (!view %in% names(VIEWS)) {
        stop(paste("unknown view:", view,
                   "(expected chimera, tecount or telocal)"))
    }
    VIEWS[[view]]
}

read_counts <- function(path) {
    m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
    storage.mode(m) <- "integer"
    m
}

transform_matrix <- function(counts, method) {
    if (method == "log2") {
        return(log2(counts + 1))
    }
    if (ncol(counts) < 2) {
        message(sprintf(
            "only %d sample(s); vst/rlog requires >= 2, falling back to log2",
            ncol(counts)
        ))
        return(log2(counts + 1))
    }
    if (method == "vst") {
        tryCatch(
            vst(counts, blind = TRUE),
            error = function(e) transform_fallback(counts, method, e)
        )
    } else if (method == "rlog") {
        tryCatch(
            rlog(counts, blind = TRUE),
            error = function(e) transform_fallback(counts, method, e)
        )
    } else {
        stop(paste("unknown transform:", method))
    }
}

transform_fallback <- function(counts, method, orig_error) {
    # vst()/rlog() fail on small matrices (fewer rows than DESeq2's nsub);
    # the direct *Transformation() calls handle those but can still choke when
    # every event's gene-wise dispersion sits at the minimum (curve fitting
    # gives up), which happens easily with a handful of chimeric junctions.
    direct <- if (method == "vst") {
        tryCatch(
            varianceStabilizingTransformation(counts, blind = TRUE),
            error = function(e) NULL
        )
    } else {
        tryCatch(
            rlogTransformation(counts, blind = TRUE),
            error = function(e) NULL
        )
    }
    if (!is.null(direct)) {
        return(as.matrix(direct))
    }
    message(sprintf(
        "%s unavailable for %d features (both vst/rlog and the direct varianceStabilizingTransformation/rlogTransformation failed: %s); falling back to log2(counts + 1) for the QC view",
        method, nrow(counts), conditionMessage(orig_error)
    ))
    log2(counts + 1)
}

load_conditions <- function(samples_path) {
    sheet <- read.csv(samples_path, comment.char = "#", check.names = FALSE)
    if ("condition" %in% colnames(sheet)) {
        cond <- sheet$condition
    } else {
        cond <- rep("all", nrow(sheet))
    }
    cond[is.na(cond) | cond == ""] <- "all"
    setNames(as.character(cond), as.character(sheet$sample))
}

# --- MultiQC custom-content JSON writers -------------------------------------
# The plots mode emits two documents MultiQC picks up by the `_mqc.json`
# suffix (results/chimera/qc/*_mqc.json) and renders as an interactive
# scatter (PCA) and heatmap (sample distances). Hand-built JSON: the payloads
# are small and fixed, and the QC env needs no extra JSON dependency.

json_escape <- function(x) gsub('([\\\\"])', '\\\\\\1', x)
json_num <- function(x) sprintf("%.6g", as.numeric(x))
json_str_arr <- function(x) paste(sprintf('"%s"', json_escape(x)), collapse = ", ")

write_pca_mqc <- function(path, samples, x, y, colors, pc1, pc2, transform,
                          v, note = NULL) {
    pts <- vapply(seq_along(samples), function(i) {
        sprintf('"%s": {"x": %s, "y": %s, "color": "%s"}',
                json_escape(samples[i]), json_num(x[i]), json_num(y[i]), colors[i])
    }, character(1))
    desc <- if (is.null(note)) {
        sprintf(paste0("Principal-component analysis of the %s-transformed ",
                       "%s counts matrix, colored by sample condition."),
                transform, v$id)
    } else {
        note
    }
    body <- paste0(
        '{\n',
        sprintf('  "id": "%s_chimera_reads_sample_qc_pca",\n', v$id),
        sprintf('  "parent_id": "%s",\n', v$parent),
        sprintf('  "parent_name": "%s",\n', v$label),
        sprintf('  "section_name": "%sPCA",\n', v$section_prefix),
        sprintf('  "description": "%s",\n', json_escape(desc)),
        '  "plot_type": "scatter",\n',
        '  "pconfig": {\n',
        sprintf('    "id": "%s_pca_plot",\n', v$id),
        sprintf('    "title": "%s counts: PCA",\n', v$label),
        sprintf('    "xlab": "PC1 (%s%%)",\n', json_num(pc1)),
        sprintf('    "ylab": "PC2 (%s%%)"\n', json_num(pc2)),
        '  },\n',
        sprintf('  "data": {%s}\n', paste(pts, collapse = ", ")),
        '}\n'
    )
    writeLines(body, path)
}

write_heatmap_mqc <- function(path, samples, d, transform, v, note = NULL) {
    desc <- if (is.null(note)) {
        sprintf(paste0("Pairwise Euclidean distances between samples on the ",
                       "%s-transformed %s counts matrix."),
                transform, v$id)
    } else {
        note
    }
    rows <- apply(d, 1, function(r) paste0('[', paste(json_num(r), collapse = ", "), ']'))
    body <- paste0(
        '{\n',
        sprintf('  "id": "%s_chimera_reads_sample_qc_heatmap",\n', v$id),
        sprintf('  "parent_id": "%s",\n', v$parent),
        sprintf('  "parent_name": "%s",\n', v$label),
        sprintf('  "section_name": "%sClusters",\n', v$section_prefix),
        sprintf('  "description": "%s",\n', json_escape(desc)),
        '  "plot_type": "heatmap",\n',
        '  "pconfig": {\n',
        sprintf('    "id": "%s_heatmap_plot",\n', v$id),
        '    "title": "Euclidean distance between samples"\n',
        '  },\n',
        sprintf('  "ycats": [%s],\n', json_str_arr(samples)),
        sprintf('  "xcats": [%s],\n', json_str_arr(samples)),
        sprintf('  "data": [%s]\n', paste(rows, collapse = ",\n    ")),
        '}\n'
    )
    writeLines(body, path)
}

write_empty_mqc <- function(path, kind, v) {
    # Valid custom-content documents for the no-data case, so the multiqc rule
    # always sees its inputs and the report still documents why nothing is
    # plotted.
    note <- sprintf("No %s passed the QC-view filters; plot skipped.", v$noun_plural)
    body <- if (kind == "scatter") {
        paste0(
            '{\n',
            sprintf('  "id": "%s_chimera_reads_sample_qc_pca",\n', v$id),
            sprintf('  "parent_id": "%s",\n', v$parent),
            sprintf('  "parent_name": "%s",\n', v$label),
            sprintf('  "section_name": "%sPCA",\n', v$section_prefix),
            sprintf('  "description": "%s",\n', json_escape(note)),
            '  "plot_type": "scatter",\n',
            sprintf('  "pconfig": {"id": "%s_pca_plot", "title": "%s counts: PCA"},\n',
                    v$id, v$label),
            '  "data": {}\n',
            '}\n'
        )
    } else {
        paste0(
            '{\n',
            sprintf('  "id": "%s_chimera_reads_sample_qc_heatmap",\n', v$id),
            sprintf('  "parent_id": "%s",\n', v$parent),
            sprintf('  "parent_name": "%s",\n', v$label),
            sprintf('  "section_name": "%sClusters",\n', v$section_prefix),
            sprintf('  "description": "%s",\n', json_escape(note)),
            '  "plot_type": "heatmap",\n',
            sprintf('  "pconfig": {"id": "%s_heatmap_plot", "title": "Euclidean distance between samples"},\n',
                    v$id),
            '  "ycats": [],\n',
            '  "xcats": [],\n',
            '  "data": []\n',
            '}\n'
        )
    }
    writeLines(body, path)
}

do_transform <- function(argv) {
    v <- view_params(argv[2])
    counts_path <- argv[3]
    samples_path <- argv[4]
    method <- argv[5]
    min_samples_present <- as.integer(argv[6])
    min_total_counts <- as.integer(argv[7])
    out_matrix <- argv[8]

    counts <- read_counts(counts_path)
    if (nrow(counts) == 0) {
        file.create(out_matrix)
        message(sprintf("no %s; wrote empty transformed matrix", v$noun_plural))
        quit(save = "no", status = 0)
    }
    # QC-view filter: keep features seen in >= min_samples_present samples with
    # >= min_total_counts supporting reads overall. Only affects this view.
    present <- rowSums(counts > 0)
    keep <- present >= min_samples_present &
        rowSums(counts) >= min_total_counts
    counts <- counts[keep, , drop = FALSE]
    if (nrow(counts) == 0) {
        file.create(out_matrix)
        message(sprintf("no %s passed the QC-view filters; empty transformed matrix",
                        v$noun_plural))
        quit(save = "no", status = 0)
    }
    message(sprintf("transform=%s on %d %s x %d samples", method,
                    nrow(counts), v$noun_plural, ncol(counts)))
    transformed <- as.matrix(transform_matrix(counts, method))
    write.table(transformed, out_matrix, sep = "\t", quote = FALSE)
}

do_plots <- function(argv) {
    v <- view_params(argv[2])
    matrix_path <- argv[3]
    samples_path <- argv[4]
    min_events <- as.integer(argv[5])
    transform <- argv[6]
    out_pca <- argv[7]
    out_heatmap <- argv[8]

    sz <- file.info(matrix_path)$size
    tmat <- if (!is.na(sz) && sz > 0) {
        as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
    } else {
        # the transform writes a 0-byte file when nothing passed its QC filter
        matrix(nrow = 0, ncol = 0)
    }
    empty <- function(msg) {
        message(msg)
        write_empty_mqc(out_pca, "scatter", v)
        write_empty_mqc(out_heatmap, "heatmap", v)
        quit(save = "no", status = 0)
    }
    if (ncol(tmat) < 2) {
        # a single sample can't support a PCA or a sample-distance view
        empty("only one sample; skipping PCA/heatmap")
    }
    # Features that transform to a constant (or NaN/Inf) column carry no signal
    # and make prcomp(scale.=TRUE) fail outright, so drop them for the QC view.
    finite_rows <- apply(tmat, 1, function(row) all(is.finite(row)))
    var_rows <- apply(tmat, 1, var, na.rm = TRUE) > 0
    dropped <- sum(!(finite_rows & var_rows))
    if (dropped > 0) {
        message(sprintf("dropping %d zero-variance/non-finite %s(s) from the QC view",
                        dropped, v$noun_singular))
    }
    tmat <- tmat[finite_rows & var_rows, , drop = FALSE]
    if (nrow(tmat) < max(min_events, 2)) {
        # too few features for a meaningful PCA/distance view; ship empty custom
        # content rather than a cryptic prcomp error.
        empty(sprintf("fewer than %d %s; skipping PCA/heatmap", min_events,
                      v$noun_plural))
    }

    conditions <- load_conditions(samples_path)
    samples <- colnames(tmat)
    groups <- conditions[samples]
    groups[is.na(groups)] <- "all"

    pca <- prcomp(t(tmat), center = TRUE, scale. = TRUE)
    pc1 <- round(100 * summary(pca)$importance[2, 1], 1)
    pc2 <- round(100 * summary(pca)$importance[2, 2], 1)

    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7")
    ugroups <- unique(groups)
    group_colors <- setNames(palette[seq_along(ugroups)], ugroups)
    point_colors <- unname(group_colors[groups])

    write_pca_mqc(out_pca, samples, pca$x[, 1], pca$x[, 2], point_colors,
                  pc1, pc2, transform, v)

    write_heatmap_mqc(out_heatmap, samples, as.matrix(dist(t(tmat))), transform, v)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("usage: sample_qc.R --transform view ... | --plots view ...")
}
mode <- args[1]
if (mode == "--transform") {
    do_transform(args)
} else if (mode == "--plots") {
    do_plots(args)
} else {
    stop(paste("unknown mode:", mode))
}
