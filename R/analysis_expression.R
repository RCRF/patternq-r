# Broad-scale descriptive expression analysis, backported from the Clojure
# variant-forensics reports (unify-central/analysis):
#
#   - a sample vs a reference cohort: z-scores, percentile ranks, top genes
#     (quantitative.clj rank-expression, util.clj percentile-rank)
#   - geneset views of samples against cohort distributions
#     (util.clj examine-geneset*, plot-genex-vs-cohort*)
#   - two-sample change: log fold change and MA plots
#     (quantitative.clj log-fold-change, H37001 gene_expression_change.clj)
#   - ssGSEA, top varying genes, expression distances / nearest samples
#
# Everything here works on data returned by gene_expression(); the only
# queries are gene_expression() calls (batched over genes for cohorts).

#' Gene sets shipped with patternq
#'
#' Hallmark sets (apoptosis, DNA repair, hypoxia, inflammatory), antibody
#' therapy targets, housekeeping genes, germline multi-cancer panel, melanoma
#' phenotype sets, neural crest and adult kidney reference sets.
#'
#' @param name Gene set name (see \code{names(genesets())})
#' @return \code{genesets()}: named list of character vectors;
#'   \code{geneset(name)}: character vector
#' @export
genesets <- function() {
    dir <- system.file("genesets", package = "patternq")
    files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
    out <- lapply(files, function(f) {
        g <- trimws(readLines(f, warn = FALSE))
        g[nzchar(g)]
    })
    names(out) <- sub("\\.txt$", "", basename(files))
    out
}

#' @rdname genesets
#' @export
geneset <- function(name) {
    gs <- genesets()
    if (!name %in% names(gs)) stop(sprintf("Unknown gene set '%s'", name), call. = FALSE)
    gs[[name]]
}

#' Expression of one sample as a named vector
#'
#' @inheritParams gene_expression
#' @param sample Sample id
#' @return Named numeric vector (HGNC symbol -> value; values of several gene
#'   products of the same gene are summed)
#' @export
sample_expression <- function(sample, db = NULL, measurement = "tpm", measurement.set = NULL, genes = NULL) {
    gx <- gene_expression(db = db, genes = genes, samples = sample, measurement = measurement,
                          measurement.set = measurement.set)
    if (!nrow(gx)) return(stats::setNames(numeric(0), character(0)))
    v <- tapply(gx$value, gx$hgnc_symbol, sum)
    keep_provenance(stats::setNames(as.numeric(v), names(v)), gx)
}

#' Percentile rank of a value within a distribution
#'
#' Mid-rank percentile: (count below + half the ties) / n, as a percentage.
#'
#' @param values Reference values
#' @param x Value(s) to rank
#' @return Percentile(s) in [0, 100]
#' @export
percentile_rank <- function(values, x) {
    values <- values[!is.na(values)]
    n <- length(values)
    if (!n) return(rep(NA_real_, length(x)))
    vapply(x, function(v) 100 * (sum(values < v) + 0.5 * sum(values == v)) / n, numeric(1))
}

#' Compare a sample's expression to a reference cohort
#'
#' For every gene, the sample's value is placed within the cohort's
#' distribution: cohort mean/sd/median on \code{log2(1 + x)} scale (or raw with
#' \code{log = FALSE}), z-score, and percentile rank. The cohort is streamed in
#' gene batches, so a whole transcriptome against a large cohort never needs to
#' be held in memory.
#'
#' Measurements must be comparable (same units / normalization) between the
#' sample's database and the cohort's; TPM is the usual common ground.
#'
#' @param sample Sample id in \code{db}
#' @param db Database holding the sample
#' @param cohort.db Reference cohort database, e.g. \code{resolve_db("tcga-uvm")}
#' @param measurement Measurement attribute for the sample
#' @param cohort.measurement Measurement attribute for the cohort (defaults to
#'   \code{measurement})
#' @param genes Genes to compare (default: all genes measured in the sample)
#' @param measurement.set,cohort.measurement.set Optional measurement sets
#' @param log Compare on \code{log2(1 + x)}
#' @param fill.missing Count a gene with no stored value in a cohort sample as
#'   0. Imports often omit zero measurements, so without this a gene's cohort
#'   distribution is built from only the samples that express it. Cohort
#'   samples are those with a value for \code{anchor.gene}.
#' @param anchor.gene Housekeeping gene that defines the cohort's samples
#' @param sd.floor Minimum cohort standard deviation (on the comparison scale)
#'   so near-constant genes do not produce enormous z-scores
#' @param batch.size Genes per cohort query
#' @param min.cohort Minimum cohort values for a gene to be scored
#' @return \code{data.frame}: \code{hgnc_symbol}, \code{value}, \code{z},
#'   \code{percentile}, \code{cohort_n} (cohort size), \code{cohort_observed}
#'   (cohort samples with a stored value; low coverage of a gene that is
#'   usually expressed points to an import or annotation problem),
#'   \code{cohort_mean}, \code{cohort_sd}, \code{cohort_median}
#' @export
compare_to_cohort <- function(sample, db = NULL, cohort.db, measurement = "tpm",
                              cohort.measurement = measurement, genes = NULL,
                              measurement.set = NULL, cohort.measurement.set = NULL,
                              log = TRUE, fill.missing = TRUE, anchor.gene = "GAPDH",
                              sd.floor = 0.25, batch.size = 2000, min.cohort = 5) {
    x <- sample_expression(sample, db = db, measurement = measurement,
                           measurement.set = measurement.set, genes = genes)
    if (!length(x)) stop(sprintf("No %s expression for sample %s", measurement, sample), call. = FALSE)
    tr <- if (log) function(v) log2(1 + pmax(v, 0)) else identity
    cohort.samples <- NULL
    if (fill.missing) {
        anchor <- gene_expression(db = cohort.db, genes = anchor.gene, measurement = cohort.measurement,
                                  measurement.set = cohort.measurement.set)
        cohort.samples <- unique(anchor$sample_id)
        if (!length(cohort.samples))
            stop(sprintf("No %s values for %s in %s; pass anchor.gene or fill.missing = FALSE",
                         cohort.measurement, anchor.gene, cohort.db), call. = FALSE)
    }
    batches <- split(names(x), ceiling(seq_along(x) / batch.size))
    parts <- lapply(batches, function(b) {
        cg <- gene_expression(db = cohort.db, genes = b, measurement = cohort.measurement,
                              measurement.set = cohort.measurement.set)
        if (!nrow(cg)) return(NULL)
        cg <- stats::aggregate(value ~ sample_id + hgnc_symbol, data = cg, FUN = sum)
        by.gene <- split(tr(cg$value), cg$hgnc_symbol)
        observed <- vapply(by.gene, length, integer(1))
        if (!is.null(cohort.samples)) {
            n.all <- length(cohort.samples)
            by.gene <- lapply(by.gene, function(v) c(v, rep(tr(0), max(0, n.all - length(v)))))
        }
        g <- names(by.gene)
        v <- tr(x[g])
        data.frame(hgnc_symbol = g, value = unname(x[g]),
                   cohort_n = vapply(by.gene, length, integer(1)),
                   cohort_observed = observed,
                   cohort_mean = vapply(by.gene, mean, numeric(1)),
                   cohort_sd = vapply(by.gene, stats::sd, numeric(1)),
                   cohort_median = vapply(by.gene, stats::median, numeric(1)),
                   percentile = mapply(function(vals, xv) percentile_rank(vals, xv), by.gene, v),
                   row.names = NULL, stringsAsFactors = FALSE)
    })
    out <- bind_rows_fill(parts)
    out <- out[out$cohort_n >= min.cohort, , drop = FALSE]
    sdv <- pmax(ifelse(is.na(out$cohort_sd), 0, out$cohort_sd), sd.floor)
    out$z <- (tr(out$value) - out$cohort_mean) / sdv
    out <- order_columns(out, c("hgnc_symbol", "value", "z", "percentile", "cohort_n", "cohort_observed", "cohort_mean",
                                "cohort_sd", "cohort_median"))
    out <- out[order(-abs(out$z)), , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "patternq_comparison") <- list(sample = sample, db = ensure_db(db), cohort.db = cohort.db,
                                             measurement = measurement, cohort.measurement = cohort.measurement,
                                             log = log, cohort.size = length(cohort.samples))
    out
}

#' Top genes by z-score against a cohort
#'
#' @param comparison Output of \code{\link{compare_to_cohort}}
#' @param n Number of genes
#' @param direction \code{"both"} (largest |z|), \code{"up"} or \code{"down"}
#' @param min.value Ignore genes whose sample value and cohort median are both
#'   below this (removes noise among barely-expressed genes)
#' @param min.observed Minimum fraction of the cohort with a stored value for
#'   the gene (guards against genes missing or mis-annotated in the cohort's
#'   import dominating the ranking)
#' @return Rows of \code{comparison}
#' @export
top_by_zscore <- function(comparison, n = 25, direction = c("both", "up", "down"), min.value = 1,
                          min.observed = 0.5) {
    direction <- match.arg(direction)
    cmp <- comparison[!is.na(comparison$z), , drop = FALSE]
    if ("cohort_observed" %in% names(cmp))
        cmp <- cmp[cmp$cohort_observed >= min.observed * cmp$cohort_n, , drop = FALSE]
    expressed <- cmp$value >= min.value | (2^cmp$cohort_median - 1) >= min.value
    cmp <- cmp[expressed, , drop = FALSE]
    o <- switch(direction, both = order(-abs(cmp$z)), up = order(-cmp$z), down = order(cmp$z))
    cmp <- cmp[o, , drop = FALSE]
    if (direction == "up") cmp <- cmp[cmp$z > 0, , drop = FALSE]
    if (direction == "down") cmp <- cmp[cmp$z < 0, , drop = FALSE]
    utils::head(cmp, n)
}

#' Plot top genes by z-score
#'
#' Horizontal diverging bars of z-scores vs the reference cohort.
#'
#' @param comparison Output of \code{\link{compare_to_cohort}} (or of
#'   \code{\link{top_by_zscore}})
#' @param n Number of genes (largest |z|)
#' @param min.value,min.observed See \code{\link{top_by_zscore}}
#' @param title Plot title
#' @return A plotly object
#' @export
plot_zscores <- function(comparison, n = 30, min.value = 1, min.observed = 0.5, title = NULL) {
    top <- top_by_zscore(comparison, n = n, min.value = min.value, min.observed = min.observed)
    top <- top[order(top$z), , drop = FALSE]
    d <- pq_theme()$diverging
    info <- attr(comparison, "patternq_comparison")
    if (is.null(title) && !is.null(info))
        title <- sprintf("%s vs %s: top genes by z-score", info$sample, info$cohort.db)
    p <- plotly::plot_ly(x = top$z, y = top$hgnc_symbol, type = "bar", orientation = "h",
                         marker = list(color = ifelse(top$z >= 0, d$high, d$low)),
                         text = sprintf("value %.3g \u00b7 cohort median %.3g \u00b7 %.1f pct",
                                        top$value, 2^top$cohort_median - 1, top$percentile),
                         textposition = "none",
                         hovertemplate = "%{y}<br>z = %{x:.2f}<br>%{text}<extra></extra>")
    pq_layout(p, title = title, showlegend = FALSE, bargap = 0.25,
              xaxis = list(title = "z-score vs cohort (log2(1+x))", zeroline = TRUE),
              yaxis = list(title = "", type = "category", categoryorder = "array",
                           categoryarray = top$hgnc_symbol))
}

#' Samples against cohort distributions for a gene set
#'
#' One row per gene: the cohort's distribution (violin or box, one per
#' cohort) with each sample's value marked. Log scale by default (values are
#' floored at \code{floor} so zeros can be drawn).
#'
#' @param sample.expr Long data from \code{\link{gene_expression}} for the
#'   samples (\code{sample_id}, \code{hgnc_symbol}, \code{value})
#' @param cohort.expr Long data for the cohort(s); a \code{cohort} column
#'   distinguishes several cohorts (see \code{\link{examine_geneset}})
#' @param type \code{"violin"} or \code{"box"}
#' @param log Log-scaled value axis
#' @param floor Minimum drawn value on log scale
#' @param title,xlab Labels
#' @return A plotly object
#' @export
plot_vs_cohort <- function(sample.expr, cohort.expr, type = c("violin", "box"), log = TRUE,
                           floor = 0.01, title = "Samples vs cohort expression", xlab = "expression") {
    type <- match.arg(type)
    if (!"cohort" %in% names(cohort.expr)) cohort.expr$cohort <- "cohort"
    cohort.expr <- stats::aggregate(value ~ cohort + sample_id + hgnc_symbol, data = cohort.expr, FUN = sum)
    sample.expr <- stats::aggregate(value ~ sample_id + hgnc_symbol, data = sample.expr, FUN = sum)
    genes <- unique(c(unique(sample.expr$hgnc_symbol), unique(cohort.expr$hgnc_symbol)))
    # on a log scale, draw log10 values on a linear axis (so violin densities
    # are estimated on the log scale) and label ticks as powers of ten
    fl <- function(v) if (log) log10(pmax(v, floor)) else v
    cohorts <- sort(unique(cohort.expr$cohort))
    ccols <- series_colors(length(cohorts))
    th <- pq_theme()
    p <- plotly::plot_ly(height = 160 + 42 * length(genes))
    for (i in seq_along(cohorts)) {
        ce <- cohort.expr[cohort.expr$cohort == cohorts[i], ]
        if (type == "violin")
            p <- plotly::add_trace(p, type = "violin", orientation = "h", x = fl(ce$value), y = ce$hgnc_symbol,
                                   name = cohorts[i], legendgroup = cohorts[i],
                                   line = list(color = ccols[i], width = 1), fillcolor = alpha_color(ccols[i], 0.25),
                                   points = FALSE, spanmode = "hard", scalemode = "width", width = 0.8,
                                   box = list(visible = TRUE, fillcolor = th$surface, line = list(color = ccols[i]), width = 0.2),
                                   meanline = list(visible = FALSE), hoverinfo = "y+name")
        else
            p <- plotly::add_boxplot(p, orientation = "h", x = fl(ce$value), y = ce$hgnc_symbol,
                                     name = cohorts[i], legendgroup = cohorts[i],
                                     marker = list(color = ccols[i], size = 3), line = list(color = ccols[i]),
                                     fillcolor = alpha_color(ccols[i], 0.2), boxpoints = "outliers")
    }
    symbols <- c("diamond", "square", "circle", "triangle-up", "x", "star")
    sids <- sort(unique(sample.expr$sample_id))
    for (j in seq_along(sids)) {
        se <- sample.expr[sample.expr$sample_id == sids[j], ]
        p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = fl(se$value), y = se$hgnc_symbol,
                               name = sids[j],
                               marker = list(symbol = symbols[(j - 1) %% length(symbols) + 1], size = 11,
                                             color = th$text_primary, line = list(color = th$surface, width = 1.5)),
                               customdata = se$value,
                               hovertemplate = paste0(sids[j], "<br>%{y}: %{customdata:.3g}<extra></extra>"))
    }
    xaxis <- list(title = xlab)
    if (log) {
        rng <- range(c(fl(cohort.expr$value), fl(sample.expr$value)), finite = TRUE)
        ticks <- seq(floor(rng[1]), ceiling(rng[2]))
        xaxis$tickvals <- ticks
        xaxis$ticktext <- ifelse(ticks >= 3, formatC(10^ticks, format = "d", big.mark = ","),
                                 format(10^ticks, scientific = FALSE, drop0trailing = TRUE))
    }
    pq_layout(p, title = title, boxmode = "group", xaxis = xaxis,
              yaxis = list(title = "", type = "category", categoryorder = "array",
                           categoryarray = rev(genes), autorange = TRUE))
}

#' Examine a gene set: samples vs reference cohort(s)
#'
#' Fetches expression for \code{genes} in the samples and in each cohort and
#' draws \code{\link{plot_vs_cohort}}. The port of the Clojure reports'
#' \code{examine-geneset} family.
#'
#' @param genes HGNC symbols (or a name from \code{\link{genesets}})
#' @param samples Sample ids in \code{db}
#' @param db Database holding the samples
#' @param cohort.dbs Reference cohort database(s); names of the vector (if
#'   any) label the cohorts
#' @param measurement,cohort.measurement Measurement attributes
#' @param type \code{"violin"} or \code{"box"}
#' @param title Plot title
#' @param ... Passed to \code{\link{plot_vs_cohort}}
#' @return A plotly object
#' @export
examine_geneset <- function(genes, samples, db = NULL, cohort.dbs, measurement = "tpm",
                            cohort.measurement = measurement, type = c("violin", "box"),
                            title = NULL, ...) {
    if (length(genes) == 1 && genes %in% names(genesets())) {
        if (is.null(title)) title <- genes
        genes <- geneset(genes)
    }
    se <- gene_expression(db = db, genes = genes, samples = samples, measurement = measurement)
    labels <- if (is.null(names(cohort.dbs))) cohort.dbs else names(cohort.dbs)
    ce <- bind_rows_fill(lapply(seq_along(cohort.dbs), function(i) {
        d <- gene_expression(db = cohort.dbs[[i]], genes = genes, measurement = cohort.measurement)
        if (nrow(d)) d$cohort <- labels[i]
        d
    }))
    plot_vs_cohort(se, ce, type = match.arg(type),
                   title = if (is.null(title)) "Samples vs cohort expression" else title,
                   xlab = paste(measurement, "(log scale)"), ...)
}

#' Log fold change between two samples
#'
#' \code{lfc = log2((b + pseudocount) / (a + pseudocount))}; \code{A} (the MA
#' plot's x) is the mean of \code{log10(1 + a)} and \code{log10(1 + b)}, as in
#' the Clojure reports.
#'
#' @param a,b Named numeric vectors (e.g. from \code{\link{sample_expression}});
#'   genes missing in one are 0
#' @param pseudocount Added before the ratio
#' @return \code{data.frame}: \code{hgnc_symbol}, \code{value_a},
#'   \code{value_b}, \code{lfc}, \code{avg_log10}
#' @export
log_fold_change <- function(a, b, pseudocount = 1) {
    genes <- union(names(a), names(b))
    va <- unname(a[genes]); va[is.na(va)] <- 0
    vb <- unname(b[genes]); vb[is.na(vb)] <- 0
    data.frame(hgnc_symbol = genes, value_a = va, value_b = vb,
               lfc = log2((vb + pseudocount) / (va + pseudocount)),
               avg_log10 = (log10(1 + va) + log10(1 + vb)) / 2,
               stringsAsFactors = FALSE)
}

#' Compare expression between two samples
#'
#' @param sample.a Reference sample (e.g. baseline)
#' @param sample.b Comparison sample (e.g. later timepoint)
#' @param db Database (both samples); use \code{db.b} if \code{sample.b}
#'   lives elsewhere
#' @param db.b Database of \code{sample.b}
#' @param measurement Measurement attribute
#' @param min.avg Drop genes with \code{avg_log10} below this (low expression)
#' @param pseudocount See \code{\link{log_fold_change}}
#' @return \code{\link{log_fold_change}} output, ordered by \code{lfc}
#' @export
compare_samples <- function(sample.a, sample.b, db = NULL, db.b = db, measurement = "tpm",
                            min.avg = 0.5, pseudocount = 1) {
    a <- sample_expression(sample.a, db = db, measurement = measurement)
    b <- sample_expression(sample.b, db = db.b, measurement = measurement)
    out <- log_fold_change(a, b, pseudocount = pseudocount)
    out <- out[out$avg_log10 >= min.avg, , drop = FALSE]
    out <- out[order(-out$lfc), , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "patternq_comparison") <- list(sample.a = sample.a, sample.b = sample.b, measurement = measurement)
    out
}

#' MA plot of expression change between two samples
#'
#' x: average expression (\code{avg_log10}); y: log2 fold change. Points past
#' \code{lfc.threshold} are colored (up red, down blue, on the diverging
#' scale), the rest recede; \code{highlight} genes and the top
#' \code{label.top} each way are labelled.
#'
#' @param change Output of \code{\link{compare_samples}} or \code{\link{log_fold_change}}
#' @param lfc.threshold Fold change (log2) guides and coloring
#' @param highlight Genes to label
#' @param label.top Label this many top up and down genes (above
#'   \code{label.min.avg})
#' @param label.min.avg Minimum \code{avg_log10} for automatic labels
#' @param title Plot title
#' @return A plotly object
#' @export
plot_ma <- function(change, lfc.threshold = 2.5, highlight = NULL, label.top = 8,
                    label.min.avg = 1.5, title = NULL) {
    d <- pq_theme()$diverging
    th <- pq_theme()
    info <- attr(change, "patternq_comparison")
    if (is.null(title) && !is.null(info$sample.a))
        title <- sprintf("Expression change: %s \u2192 %s", info$sample.a, info$sample.b)
    cls <- ifelse(change$lfc >= lfc.threshold, "up", ifelse(change$lfc <= -lfc.threshold, "down", "unchanged"))
    p <- plotly::plot_ly()
    spec <- list(unchanged = list(col = alpha_color("#8a8983", 0.35), name = "within threshold"),
                 down = list(col = d$low, name = sprintf("down (lfc \u2264 -%s)", lfc.threshold)),
                 up = list(col = d$high, name = sprintf("up (lfc \u2265 %s)", lfc.threshold)))
    for (k in c("unchanged", "down", "up")) {
        s <- change[cls == k, , drop = FALSE]
        if (!nrow(s)) next
        p <- plotly::add_trace(p, type = "scattergl", mode = "markers", x = s$avg_log10, y = s$lfc,
                               text = s$hgnc_symbol, name = spec[[k]]$name,
                               marker = list(color = spec[[k]]$col, size = if (k == "unchanged") 5 else 7),
                               customdata = cbind(s$value_a, s$value_b),
                               hovertemplate = "%{text}<br>lfc %{y:.2f}<br>a %{customdata[0]:.3g} \u2192 b %{customdata[1]:.3g}<extra></extra>")
    }
    lab <- change[change$avg_log10 >= label.min.avg, , drop = FALSE]
    lab <- rbind(utils::head(lab[order(-lab$lfc), ], label.top), utils::head(lab[order(lab$lfc), ], label.top))
    lab <- lab[abs(lab$lfc) >= lfc.threshold, , drop = FALSE]
    if (!is.null(highlight)) lab <- unique(rbind(change[change$hgnc_symbol %in% highlight, , drop = FALSE], lab))
    # greedy thinning so labels do not pile up (highlighted genes win)
    keep <- logical(nrow(lab))
    xr <- diff(range(change$avg_log10)); yr <- diff(range(change$lfc))
    for (i in seq_len(nrow(lab))) {
        near <- which(keep)[abs(lab$avg_log10[keep] - lab$avg_log10[i]) < 0.06 * xr &
                            abs(lab$lfc[keep] - lab$lfc[i]) < 0.05 * yr]
        keep[i] <- length(near) == 0
    }
    lab <- lab[keep, , drop = FALSE]
    ann <- lapply(seq_len(nrow(lab)), function(i) list(x = lab$avg_log10[i], y = lab$lfc[i], text = lab$hgnc_symbol[i],
                                                       showarrow = TRUE, arrowhead = 0, arrowwidth = 1,
                                                       arrowcolor = th$text_secondary, ax = 18,
                                                       ay = if (lab$lfc[i] > 0) -16 else 16,
                                                       font = list(size = 11, color = th$text_primary)))
    guide <- function(y) list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = y, y1 = y,
                              line = list(color = th$text_secondary, width = 1, dash = if (y == 0) "solid" else "dot"))
    pq_layout(p, title = title, annotations = ann,
              shapes = list(guide(0), guide(lfc.threshold), guide(-lfc.threshold)),
              hovermode = "closest",
              xaxis = list(title = "average expression, log10(1 + x)"),
              yaxis = list(title = "log2 fold change"))
}

#' Fold change bars for selected genes
#'
#' @param change Output of \code{\link{compare_samples}}
#' @param genes Genes to show (default: top 15 up and 15 down)
#' @param title Plot title
#' @return A plotly object
#' @export
plot_fold_change <- function(change, genes = NULL, title = "Change in gene expression") {
    if (is.null(genes)) {
        s <- change[order(change$lfc), ]
        genes <- unique(c(utils::head(s$hgnc_symbol, 15), utils::tail(s$hgnc_symbol, 15)))
    }
    s <- change[change$hgnc_symbol %in% genes, , drop = FALSE]
    s <- s[order(s$lfc), , drop = FALSE]
    d <- pq_theme()$diverging
    p <- plotly::plot_ly(x = s$lfc, y = s$hgnc_symbol, type = "bar", orientation = "h", height = 140 + 22 * nrow(s),
                         marker = list(color = ifelse(s$lfc >= 0, d$high, d$low)),
                         hovertemplate = "%{y}: %{x:.2f}<extra></extra>")
    pq_layout(p, title = title, showlegend = FALSE,
              xaxis = list(title = "log2 fold change"),
              yaxis = list(title = "", type = "category", categoryorder = "array", categoryarray = s$hgnc_symbol))
}

#' Single-sample GSEA enrichment score
#'
#' ssGSEA (Barbie et al. 2009): genes are ranked by expression; the score
#' sums, over the ranked list, the difference between the weighted empirical
#' distribution of the gene set (weights \code{rank^alpha}, highest expression
#' = highest rank) and the distribution of the other genes.
#'
#' Note: the Clojure report version weights by the inverse rank (highest
#' expression = rank 1); this follows the published method.
#'
#' @param x Named numeric vector (one sample)
#' @param gene.set Character vector of HGNC symbols
#' @param alpha Weight exponent
#' @return Score (numeric)
#' @export
ssgsea_score <- function(x, gene.set, alpha = 0.25) {
    x <- x[!is.na(x)]
    o <- order(x, decreasing = TRUE)
    g <- names(x)[o]
    n <- length(g)
    r <- rank(x, ties.method = "first")[o]  # highest expression -> rank n
    hit <- g %in% gene.set
    nh <- sum(hit)
    if (nh == 0 || nh == n) return(NA_real_)
    w <- ifelse(hit, r^alpha, 0)
    p.hit <- cumsum(w) / sum(w)
    p.miss <- cumsum(!hit) / (n - nh)
    sum(p.hit - p.miss)
}

#' ssGSEA scores for several samples and gene sets
#'
#' @param m Matrix, samples x genes (e.g. \code{to_matrix(gene_expression(...), col = "hgnc_symbol")})
#' @param gene.sets Named list of gene sets (default: \code{\link{genesets}()})
#' @param alpha Weight exponent
#' @return Matrix, samples x gene sets
#' @export
ssgsea <- function(m, gene.sets = genesets(), alpha = 0.25) {
    out <- sapply(gene.sets, function(gs) apply(m, 1, function(row) ssgsea_score(row, gs, alpha)))
    if (!is.matrix(out)) out <- matrix(out, nrow = nrow(m), dimnames = list(rownames(m), names(gene.sets)))
    out
}

#' Most variable genes
#'
#' Variance of \code{log2(1 + x)} across samples.
#'
#' @param m Matrix, samples x genes
#' @param n Number of genes
#' @return \code{data.frame}: \code{hgnc_symbol}, \code{mean}, \code{variance}
#' @export
top_varying_genes <- function(m, n = 500) {
    l <- log2(1 + pmax(m, 0))
    v <- apply(l, 2, stats::var, na.rm = TRUE)
    mu <- colMeans(l, na.rm = TRUE)
    o <- order(-v)
    utils::head(data.frame(hgnc_symbol = colnames(m)[o], mean = mu[o], variance = v[o],
                           row.names = NULL, stringsAsFactors = FALSE), n)
}

#' Expression distance between profiles
#'
#' @param a,b Named numeric vectors; genes missing from \code{b} count as 0
#'   (as in the Clojure reports)
#' @param method \code{"cosine"} (1 - cosine similarity) or \code{"euclidean"}
#' @return Distance
#' @export
expression_distance <- function(a, b, method = c("cosine", "euclidean")) {
    method <- match.arg(method)
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    bb <- b[names(a)]; bb[is.na(bb)] <- 0
    if (method == "euclidean") return(sqrt(sum((a - bb)^2)))
    ma <- sqrt(sum(a^2)); mb <- sqrt(sum(b^2))
    if (ma == 0 || mb == 0) return(1)
    1 - sum(a * bb) / (ma * mb)
}

#' Nearest samples to a profile
#'
#' Caution (from the Clojure reports): raw-expression nearest neighbours are
#' sensitive to batch, vendor and pipeline effects; compare within a
#' consistently processed cohort.
#'
#' @param x Named numeric vector (target profile)
#' @param m Matrix, samples x genes (cohort)
#' @param method See \code{\link{expression_distance}}
#' @param n Number of samples
#' @return \code{data.frame}: \code{sample_id}, \code{distance}
#' @export
nearest_samples <- function(x, m, method = c("cosine", "euclidean"), n = 10) {
    method <- match.arg(method)
    d <- apply(m, 1, function(row) {
        row[is.na(row)] <- 0
        expression_distance(x, row, method)
    })
    o <- order(d)
    utils::head(data.frame(sample_id = rownames(m)[o], distance = d[o], row.names = NULL,
                           stringsAsFactors = FALSE), n)
}
