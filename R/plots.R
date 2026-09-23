# plotly plots, hand-rolled (no ggplot or Bioconductor plotting). Each
# returns a plotly htmlwidget; plotly::plotly_build(p)$x has the figure data.
# The theme (theme/plotly-theme.json, shared by all three libraries) uses a
# colorblind-validated categorical order: colors follow entities in fixed
# order and are never cycled; series past 8 fold into "Other".

pq_theme <- local({
    theme <- NULL
    function() {
        if (is.null(theme)) {
            path <- system.file("plotly-theme.json", package = "patternq")
            theme <<- jsonlite::fromJSON(path)
        }
        theme
    }
})

#' The patternq plot theme
#'
#' Colors and fonts used by patternq plots (shared across the R, Python and
#' Clojure libraries).
#'
#' @return List with \code{categorical}, \code{sequential}, \code{diverging},
#'   \code{status}, \code{font}, \code{surface}, ...
#' @export
plot_theme <- function() pq_theme()

series_colors <- function(n) {
    cols <- pq_theme()$categorical
    if (n > length(cols)) cols <- c(cols, rep("#8a8983", n - length(cols)))  # beyond 8: neutral "Other"
    cols[seq_len(n)]
}

# Keep at most `max` categories, folding the rest into "Other".
fold_other <- function(x, max = 8) {
    tab <- sort(table(x), decreasing = TRUE)
    if (length(tab) <= max) return(x)
    keep <- names(tab)[seq_len(max - 1)]
    ifelse(x %in% keep, as.character(x), "Other")
}

alpha_color <- function(hex, alpha) {
    rgb <- grDevices::col2rgb(hex)
    sprintf("rgba(%d,%d,%d,%.2f)", rgb[1], rgb[2], rgb[3], alpha)
}

sequential_scale <- function() {
    s <- pq_theme()$sequential
    lapply(seq_along(s), function(i) list((i - 1) / (length(s) - 1), s[i]))
}

diverging_scale <- function() {
    d <- pq_theme()$diverging
    list(list(0, d$low), list(0.5, d$mid), list(1, d$high))
}

pq_layout <- function(p, title = NULL, ...) {
    th <- pq_theme()
    axis <- list(gridcolor = th$grid, zerolinecolor = th$grid, linecolor = th$grid,
                 tickfont = list(color = th$text_secondary))
    args <- list(...)
    for (ax in c("xaxis", "yaxis")) args[[ax]] <- utils::modifyList(axis, if (is.null(args[[ax]])) list() else args[[ax]])
    do.call(plotly::layout, c(list(p, title = list(text = title, x = 0, xanchor = "left"),
                                   font = list(family = th$font, color = th$text_primary),
                                   colorway = th$categorical,
                                   paper_bgcolor = th$surface, plot_bgcolor = th$surface,
                                   hoverlabel = list(font = list(family = th$font))),
                              args))
}

#' VAF histogram
#'
#' Histogram of variant allele frequencies for one or more samples, overlaid.
#'
#' @param variants A data frame from \code{\link{variants}} (needs
#'   \code{sample_id} and \code{vaf})
#' @param samples Optional subset of sample ids to plot
#' @param bin.size Bin width (default 0.05)
#' @param title Plot title
#' @return A plotly object
#' @export
plot_vaf_histogram <- function(variants, samples = NULL, bin.size = 0.05,
                               title = "VAF histogram") {
    if (!is.null(samples)) variants <- variants[variants$sample_id %in% samples, , drop = FALSE]
    ids <- sort(unique(variants$sample_id))
    cols <- series_colors(length(ids))
    p <- plotly::plot_ly()
    for (i in seq_along(ids)) {
        p <- plotly::add_histogram(p, x = variants$vaf[variants$sample_id == ids[i]], name = ids[i],
                                   xbins = list(start = 0, end = 1, size = bin.size),
                                   marker = list(color = cols[i], line = list(color = pq_theme()$surface, width = 1)),
                                   opacity = 0.7)
    }
    pq_layout(p, title = title, barmode = "overlay", showlegend = length(ids) > 1,
              xaxis = list(title = "VAF", range = c(0, 1)),
              yaxis = list(title = "variants"))
}

#' Gene expression bar plot
#'
#' @param expr A data frame from \code{\link{gene_expression}} (needs
#'   \code{sample_id}, \code{hgnc_symbol}, \code{value})
#' @param title Plot title
#' @param ylab Y axis label
#' @param log Use a log scale y axis
#' @return A plotly object
#' @export
plot_gene_expression <- function(expr, title = "Gene expression", ylab = "value", log = FALSE) {
    agg <- stats::aggregate(value ~ sample_id + hgnc_symbol, data = expr, FUN = sum)
    ids <- sort(unique(agg$sample_id))
    cols <- series_colors(length(ids))
    p <- plotly::plot_ly()
    for (i in seq_along(ids)) {
        a <- agg[agg$sample_id == ids[i], ]
        p <- plotly::add_bars(p, x = a$hgnc_symbol, y = a$value, name = ids[i],
                              marker = list(color = cols[i]))
    }
    pq_layout(p, title = title, barmode = "group", bargap = 0.2, bargroupgap = 0.05,
              showlegend = length(ids) > 1,
              yaxis = list(title = ylab, type = if (log) "log" else "linear"),
              xaxis = list(title = ""))
}

#' Sample overview: which subjects were measured by which measurement sets
#'
#' Heatmap of sample counts per subject (rows) and measurement set (columns).
#'
#' @param sample.assays Output of \code{\link{sample_assays}}
#' @param title Plot title
#' @return A plotly object
#' @export
plot_sample_overview <- function(sample.assays, title = "Samples per subject and measurement set") {
    sa <- unique(sample.assays[, c("subject_id", "sample_id", "measurement_set_name")])
    sa$n <- 1
    m <- to_matrix(sa, row = "subject_id", col = "measurement_set_name", value = "n", fun.aggregate = sum)
    # no samples: blank cell (surface), not the lightest ramp step
    p <- plotly::plot_ly(x = colnames(m), y = rownames(m), z = m, type = "heatmap", zmin = 0,
                         colorscale = sequential_scale(), xgap = 1, ygap = 1,
                         colorbar = list(title = "samples"),
                         hovertemplate = "%{y}<br>%{x}<br>%{z} samples<extra></extra>")
    pq_layout(p, title = title, xaxis = list(title = "", tickangle = -30),
              yaxis = list(title = "subject", autorange = "reversed"))
}

#' Values by timepoint
#'
#' Box plot of values per timepoint (ordered by relative order), optionally
#' split by a grouping column (e.g. response), with per-subject trajectories.
#'
#' @param tab Long data with \code{value}, \code{timepoint_id} and
#'   (for ordering) \code{timepoint_relative_order}, e.g. from
#'   \code{\link{add_sample_context}}
#' @param group Optional grouping column (e.g. \code{"bor"})
#' @param lines Draw per-subject trajectories (needs \code{subject_id}),
#'   colored by \code{group} when given
#' @param value Column to plot (e.g. \code{"change"} from
#'   \code{\link{change_from_baseline}})
#' @param timepoints Optional timepoints to show, in order
#' @param levels Optional group order (and so color order)
#' @param title,ylab Labels
#' @return A plotly object
#' @export
plot_by_timepoint <- function(tab, group = NULL, lines = FALSE, title = NULL, ylab = "value",
                              value = "value", timepoints = NULL, levels = NULL) {
    tab$value <- tab[[value]]
    tab <- tab[!is.na(tab$value) & !is.na(tab$timepoint_id), , drop = FALSE]
    if (!is.null(timepoints)) tab <- tab[tab$timepoint_id %in% timepoints, , drop = FALSE]
    ord <- if (!is.null(timepoints)) intersect(timepoints, tab$timepoint_id)
           else if ("timepoint_relative_order" %in% names(tab))
               unique(tab$timepoint_id[order(tab$timepoint_relative_order)]) else sort(unique(tab$timepoint_id))
    tab$timepoint_id <- factor(tab$timepoint_id, levels = ord)
    p <- plotly::plot_ly()
    if (is.null(group)) {
        if (lines && "subject_id" %in% names(tab)) {
            for (sid in unique(tab$subject_id)) {
                s <- tab[tab$subject_id == sid, ]
                s <- s[order(s$timepoint_id), ]
                p <- plotly::add_trace(p, x = s$timepoint_id, y = s$value, type = "scatter", mode = "lines",
                                       line = list(color = "rgba(82,81,78,0.25)", width = 1),
                                       hoverinfo = "text", text = sid, showlegend = FALSE)
            }
        }
        p <- plotly::add_boxplot(p, x = tab$timepoint_id, y = tab$value, name = ylab,
                                 marker = list(color = series_colors(1)), line = list(color = series_colors(1)),
                                 boxpoints = "all", jitter = 0.3, pointpos = 0, showlegend = FALSE)
    } else {
        g <- fold_other(as.character(tab[[group]]))
        lv <- if (is.null(levels)) sort(unique(g[!is.na(g)])) else intersect(levels, unique(g))
        cols <- series_colors(length(lv))
        if (lines && "subject_id" %in% names(tab)) {
            # per-subject trajectories, colored by group, drawn under the boxes
            for (i in seq_along(lv)) {
                sub <- tab[!is.na(g) & g == lv[i], ]
                for (sid in unique(sub$subject_id)) {
                    s1 <- sub[sub$subject_id == sid, ]
                    s1 <- s1[order(s1$timepoint_id), ]
                    if (nrow(s1) < 2) next
                    p <- plotly::add_trace(p, x = s1$timepoint_id, y = s1$value, type = "scatter", mode = "lines",
                                           line = list(color = alpha_color(cols[i], 0.35), width = 1),
                                           legendgroup = lv[i], showlegend = FALSE, hoverinfo = "text",
                                           text = paste(sid, lv[i], sep = "<br>"))
                }
            }
        }
        for (i in seq_along(lv)) {
            s <- tab[!is.na(g) & g == lv[i], ]
            p <- plotly::add_boxplot(p, x = s$timepoint_id, y = s$value, name = lv[i], legendgroup = lv[i],
                                     marker = list(color = cols[i], size = 4), line = list(color = cols[i]),
                                     fillcolor = alpha_color(cols[i], 0.15),
                                     boxpoints = if (lines) FALSE else "all", jitter = 0.3, pointpos = 0)
        }
        p <- plotly::layout(p, boxmode = if (lines) "overlay" else "group")
    }
    pq_layout(p, title = title, xaxis = list(title = "timepoint", type = "category",
                                             categoryorder = "array", categoryarray = ord),
              yaxis = list(title = ylab))
}

#' Values by group
#'
#' Box (or violin) plot of values per group, e.g. a measurement by best overall
#' response.
#'
#' @param tab Data with \code{value} and the group column
#' @param group Grouping column
#' @param violin Draw violins instead of boxes
#' @param title,ylab Labels
#' @return A plotly object
#' @export
plot_by_group <- function(tab, group, violin = FALSE, title = NULL, ylab = "value") {
    tab <- tab[!is.na(tab$value) & !is.na(tab[[group]]), , drop = FALSE]
    g <- fold_other(as.character(tab[[group]]))
    lv <- sort(unique(g))
    cols <- series_colors(length(lv))
    p <- plotly::plot_ly()
    for (i in seq_along(lv)) {
        y <- tab$value[g == lv[i]]
        if (violin)
            p <- plotly::add_trace(p, type = "violin", x = rep(lv[i], length(y)), y = y, name = lv[i],
                                   line = list(color = cols[i], width = 1.5), fillcolor = alpha_color(cols[i], 0.25),
                                   box = list(visible = TRUE, fillcolor = pq_theme()$surface,
                                              line = list(color = cols[i]), width = 0.15),
                                   meanline = list(visible = FALSE),
                                   points = "all", jitter = 0.4, pointpos = 0,
                                   marker = list(color = cols[i], size = 5, opacity = 0.8))
        else
            p <- plotly::add_boxplot(p, x = rep(lv[i], length(y)), y = y, name = lv[i],
                                     marker = list(color = cols[i]), line = list(color = cols[i]),
                                     boxpoints = "all", jitter = 0.3, pointpos = 0)
    }
    pq_layout(p, title = title, showlegend = FALSE, xaxis = list(title = group), yaxis = list(title = ylab))
}

#' Kaplan-Meier survival curves
#'
#' Hand-rolled Kaplan-Meier estimate (no survival package), optionally by
#' group, with censoring ticks.
#'
#' @param tab Data with time and event columns (e.g. \code{\link{subject_outcomes}})
#' @param time Time column (e.g. \code{"os"}, \code{"pfs"})
#' @param event Event column, logical or 0/1 (e.g. \code{"os_event"})
#' @param group Optional grouping column
#' @param title,xlab Labels
#' @param pvalue Annotate the log-rank p-value (when there are 2+ groups)
#' @param levels Optional group order (and so color order)
#' @return A plotly object
#' @export
plot_survival <- function(tab, time = "os", event = "os_event", group = NULL,
                          title = "Survival", xlab = time, pvalue = TRUE, levels = NULL) {
    tab <- tab[!is.na(tab[[time]]) & !is.na(tab[[event]]), , drop = FALSE]
    g <- if (is.null(group)) rep("all", nrow(tab)) else fold_other(as.character(tab[[group]]))
    keep <- !is.na(g)
    tab <- tab[keep, , drop = FALSE]; g <- g[keep]
    lv <- if (is.null(levels)) sort(unique(g)) else intersect(levels, unique(g))
    cols <- series_colors(length(lv))
    p <- plotly::plot_ly()
    for (i in seq_along(lv)) {
        km <- kaplan_meier(tab[[time]][g == lv[i]], as.logical(tab[[event]][g == lv[i]]))
        nm <- sprintf("%s (n=%d)", lv[i], sum(g == lv[i]))
        p <- plotly::add_trace(p, x = c(0, km$time), y = c(1, km$surv), type = "scatter", mode = "lines",
                               line = list(shape = "hv", color = cols[i], width = 2), name = nm,
                               hovertemplate = paste0(nm, "<br>t=%{x:.1f}<br>S=%{y:.2f}<extra></extra>"))
        cens <- km[km$n.censor > 0, ]
        if (nrow(cens))
            p <- plotly::add_trace(p, x = cens$time, y = cens$surv, type = "scatter", mode = "markers",
                                   marker = list(symbol = "line-ns-open", size = 9, color = cols[i]),
                                   name = paste(nm, "censored"), showlegend = FALSE, hoverinfo = "skip")
    }
    ann <- list()
    if (pvalue && length(lv) > 1) {
        lr <- logrank_test(tab[[time]], tab[[event]], g)
        ann <- list(list(xref = "paper", yref = "paper", x = 0.02, y = 0.04, xanchor = "left", showarrow = FALSE,
                         text = sprintf("log-rank p = %s", format.pval(lr$p, digits = 2)),
                         font = list(color = pq_theme()$text_secondary)))
    }
    pq_layout(p, title = title, showlegend = length(lv) > 1, annotations = ann,
              xaxis = list(title = xlab, rangemode = "tozero"),
              yaxis = list(title = "survival probability", range = c(0, 1.02)))
}

#' Heatmap with optional clustering
#'
#' Hand-rolled clustered heatmap (\code{stats::hclust}, no ComplexHeatmap /
#' pheatmap). With \code{scale = "row"} or \code{"column"} values are z-scored
#' and drawn on a diverging scale; otherwise sequential.
#'
#' @param m Numeric matrix
#' @param scale \code{"none"}, \code{"row"} or \code{"column"}
#' @param cluster.rows,cluster.cols Reorder by hierarchical clustering
#' @param method Linkage method for \code{hclust}
#' @param title Plot title
#' @param zlab Color bar title
#' @param col.groups Optional named vector (names = column names, values =
#'   group labels, e.g. survival status) drawn as an annotation strip above
#'   the heatmap
#' @return A plotly object
#' @export
plot_heatmap <- function(m, scale = c("none", "row", "column"), cluster.rows = TRUE,
                         cluster.cols = TRUE, method = "average", title = NULL, zlab = NULL,
                         col.groups = NULL) {
    scale <- match.arg(scale)
    m <- as.matrix(m)
    if (!nrow(m) || !ncol(m)) stop("plot_heatmap: empty matrix", call. = FALSE)
    if (scale == "row") m <- scale_rows(m)
    if (scale == "column") m <- t(scale_rows(t(m)))
    reorder <- function(x) {
        if (nrow(x) < 3) return(seq_len(nrow(x)))
        x[!is.finite(x)] <- 0
        d <- stats::dist(x)
        stats::hclust(d, method = method)$order
    }
    ro <- if (cluster.rows) reorder(m) else seq_len(nrow(m))
    co <- if (cluster.cols) reorder(t(m)) else seq_len(ncol(m))
    m <- m[ro, co, drop = FALSE]
    div <- scale != "none"
    lim <- if (div) max(abs(m), na.rm = TRUE) else NULL
    p <- plotly::plot_ly(x = colnames(m), y = rownames(m), z = m, type = "heatmap",
                         colorscale = if (div) diverging_scale() else sequential_scale(),
                         zmin = if (div) -lim else NULL, zmax = if (div) lim else NULL,
                         colorbar = list(title = if (is.null(zlab)) (if (div) "z-score" else "value") else zlab),
                         hovertemplate = "%{y}<br>%{x}<br>%{z:.3g}<extra></extra>")
    if (is.null(col.groups))
        return(pq_layout(p, title = title, xaxis = list(title = "", tickangle = -45, type = "category"),
                         yaxis = list(title = "", type = "category", autorange = "reversed")))
    grp <- unname(col.groups[colnames(m)])
    lv <- sort(unique(grp[!is.na(grp)]))
    cols <- series_colors(length(lv))
    strip <- plotly::plot_ly(x = colnames(m), y = "group", z = matrix(match(grp, lv), nrow = 1),
                             type = "heatmap", showscale = FALSE, xgap = 1,
                             colorscale = lapply(seq_along(cols), function(i)
                                 list(if (length(cols) == 1) 0 else (i - 1) / (length(cols) - 1), cols[i])),
                             zmin = 1, zmax = max(1, length(lv)), text = matrix(grp, nrow = 1),
                             hovertemplate = "%{x}<br>%{text}<extra></extra>")
    # the strip's legend: one colored label per group, above the plot
    ann <- lapply(seq_along(lv), function(i)
        list(xref = "paper", yref = "paper", x = 1, y = 1.02 + 0.045 * (length(lv) - i), xanchor = "right",
             yanchor = "bottom", showarrow = FALSE, text = paste("\u25A0", lv[i]),
             font = list(color = cols[i], size = 12)))
    out <- plotly::subplot(strip, p, nrows = 2, heights = c(0.05, 0.95), shareX = TRUE, margin = 0.005)
    pq_layout(out, title = title, annotations = ann, margin = list(t = 60 + 18 * length(lv)),
              xaxis = list(title = "", tickangle = -45, type = "category"),
              yaxis = list(title = "", showticklabels = FALSE),
              yaxis2 = list(title = "", type = "category", autorange = "reversed"))
}

scale_rows <- function(m) {
    mu <- rowMeans(m, na.rm = TRUE)
    s <- apply(m, 1, stats::sd, na.rm = TRUE)
    s[is.na(s) | s == 0] <- 1
    (m - mu) / s
}

#' Mutation landscape
#'
#' Genes (rows, most frequently mutated first) by samples (columns); cells
#' show the most severe variant impact in that gene and sample, on an ordered
#' scale (modifier < low < moderate < high).
#'
#' @param variants Output of \code{\link{variants}}
#' @param n.genes Number of genes to show
#' @param genes Optional explicit genes (overrides \code{n.genes})
#' @param title Plot title
#' @return A plotly object
#' @export
plot_mutation_landscape <- function(variants, n.genes = 25, genes = NULL, title = "Mutation landscape") {
    v <- variants[!is.na(variants$hgnc_symbol), , drop = FALSE]
    levels <- c("modifier", "low", "moderate", "high")
    has.impact <- "impact" %in% names(v) && any(!is.na(v$impact))
    # without impact annotation cells are simply mutated / not mutated
    v$severity <- if (has.impact) match(v$impact, levels) else 1
    v$severity[is.na(v$severity)] <- 1
    freq <- tapply(v$sample_id, v$hgnc_symbol, function(x) length(unique(x)))
    if (is.null(genes)) genes <- names(sort(freq, decreasing = TRUE))[seq_len(min(n.genes, length(freq)))]
    v <- v[v$hgnc_symbol %in% genes, , drop = FALSE]
    m <- to_matrix(v, row = "hgnc_symbol", col = "sample_id", value = "severity", fun.aggregate = max)
    m <- m[genes[genes %in% rownames(m)], , drop = FALSE]
    # samples ordered by mutation pattern, most frequent gene first
    pres <- !is.na(m)
    so <- do.call(order, lapply(seq_len(nrow(pres)), function(i) -pres[i, ]))
    m <- m[, so, drop = FALSE]
    s <- pq_theme()$sequential
    lab <- sprintf("%s (%d)", rownames(m), freq[rownames(m)])
    if (has.impact) {
        cs <- list(list(0, s[2]), list(0.33, s[3]), list(0.34, s[4]), list(0.66, s[5]), list(0.67, s[6]), list(1, s[7]))
        txt <- matrix(levels[m], nrow = nrow(m))
        p <- plotly::plot_ly(x = colnames(m), y = lab, z = m, type = "heatmap", text = txt,
                             colorscale = cs, zmin = 1, zmax = 4, xgap = 1, ygap = 1,
                             colorbar = list(title = "impact", tickvals = 1:4, ticktext = levels),
                             hovertemplate = "%{y}<br>%{x}<br>%{text}<extra></extra>")
    } else {
        p <- plotly::plot_ly(x = colnames(m), y = lab, z = m, type = "heatmap",
                             colorscale = list(list(0, s[5]), list(1, s[5])), showscale = FALSE,
                             xgap = 1, ygap = 1,
                             hovertemplate = "%{y}<br>%{x}<br>mutated<extra></extra>")
    }
    pq_layout(p, title = title,
              xaxis = list(title = sprintf("samples (%d)", ncol(m)), showticklabels = ncol(m) <= 40, type = "category"),
              yaxis = list(title = "", autorange = "reversed", type = "category"))
}
