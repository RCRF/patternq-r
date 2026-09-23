# Reshaping measurement tables and joining context (samples, subjects,
# variants, ...), ported from wick.

#' Long measurements to a matrix
#'
#' @param df Long data frame (e.g. from \code{\link{measurements}} or
#'   \code{\link{gene_expression}})
#' @param row Row identifier column (default \code{"sample_id"})
#' @param col Column identifier column(s); several are pasted with \code{"|"}
#' @param value Value column
#' @param fun.aggregate Aggregation for duplicate cells
#' @return Numeric matrix, rows x cols, \code{NA} where missing
#' @export
to_matrix <- function(df, col, row = "sample_id", value = "value", fun.aggregate = mean) {
    cols <- if (length(col) > 1) do.call(paste, c(df[col], sep = "|")) else df[[col]]
    rows <- df[[row]]
    vals <- df[[value]]
    agg <- stats::aggregate(vals, by = list(r = rows, c = cols), FUN = fun.aggregate)
    rn <- sort(unique(agg$r)); cn <- sort(unique(agg$c))
    m <- matrix(NA_real_, nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn))
    m[cbind(match(agg$r, rn), match(agg$c, cn))] <- agg$x
    m
}

#' Matrix to long format
#'
#' @param m Matrix (rows = samples)
#' @param row.name,col.name,value.name Output column names
#' @return Long data frame without \code{NA} cells
#' @export
to_long <- function(m, row.name = "sample_id", col.name = "target", value.name = "value") {
    df <- data.frame(rep(rownames(m), times = ncol(m)), rep(colnames(m), each = nrow(m)),
                     as.vector(m), stringsAsFactors = FALSE)
    names(df) <- c(row.name, col.name, value.name)
    df[!is.na(df[[value.name]]), , drop = FALSE]
}

#' Split measurements by measurement set
#'
#' wick's \code{group_by_assay_meas_set}: one element per measurement set,
#' optionally cast to a samples x target matrix.
#'
#' @param df Long data frame with \code{measurement_set}
#' @param col Target column(s) for the matrix
#' @param wide Cast to matrices
#' @param fun.aggregate Aggregation for duplicate cells
#' @return Named list by measurement set
#' @export
split_by_measurement_set <- function(df, col = NULL, wide = TRUE, fun.aggregate = mean) {
    parts <- split(df, df$measurement_set)
    if (!wide) return(parts)
    if (is.null(col)) col <- setdiff(names(df), c("sample_id", "measurement_set", "value"))
    lapply(parts, to_matrix, col = col, fun.aggregate = fun.aggregate)
}

#' Keep or drop targets (columns) of a measurement matrix
#' @param m Matrix
#' @param include,exclude Column names
#' @return Matrix
#' @export
select_targets <- function(m, include = NULL, exclude = NULL) {
    if (!is.null(include)) m <- m[, colnames(m) %in% include, drop = FALSE]
    if (!is.null(exclude)) m <- m[, !colnames(m) %in% exclude, drop = FALSE]
    m
}

ensure_long <- function(tab, col.name = "target") {
    if (is.matrix(tab)) return(keep_provenance(to_long(tab, col.name = col.name), tab))
    tab
}

#' Join context to a measurement table
#'
#' \code{add_sample_context} joins sample attributes (and, optionally, subject
#' attributes and outcomes) by \code{sample_id}; \code{add_subject_context}
#' joins subject attributes (and optionally outcomes) by \code{subject_id};
#' \code{add_variant_context} joins variant annotations by
#' \code{variant_id}; \code{add_cnv_context} joins CNV segments by
#' \code{cnv_id}. Matrices are converted to long format first.
#'
#' @param tab Measurement table (long data frame or samples x targets matrix)
#' @param db Database name
#' @param include.subjects Also join subject context
#' @param include.outcomes Also join \code{\link{subject_outcomes}}
#' @return Data frame
#' @name add_context
NULL

#' @rdname add_context
#' @export
add_sample_context <- function(tab, db = NULL, include.subjects = TRUE, include.outcomes = FALSE) {
    tab <- ensure_long(tab)
    smp <- samples(db)
    if (any(setdiff(names(smp), "sample_id") %in% names(tab)))
        smp <- smp[, c("sample_id", setdiff(names(smp), names(tab))), drop = FALSE]
    out <- merge(tab, smp, by = "sample_id", all.x = TRUE)
    if (include.subjects && "subject_id" %in% names(out))
        out <- add_subject_context(out, db = db, include.outcomes = include.outcomes)
    keep_provenance(out, tab)
}

#' @rdname add_context
#' @export
add_subject_context <- function(tab, db = NULL, include.outcomes = FALSE) {
    if (!"subject_id" %in% names(tab)) stop("tab has no subject_id column", call. = FALSE)
    sub <- subjects(db)
    sub <- sub[, c("subject_id", setdiff(names(sub), names(tab))), drop = FALSE]
    out <- merge(tab, sub, by = "subject_id", all.x = TRUE)
    if (include.outcomes) {
        oc <- subject_outcomes(db)
        oc <- oc[, c("subject_id", setdiff(names(oc), names(out))), drop = FALSE]
        out <- merge(out, oc, by = "subject_id", all.x = TRUE)
    }
    keep_provenance(out, tab)
}

#' @rdname add_context
#' @export
add_variant_context <- function(tab, db = NULL) {
    tab <- ensure_long(tab, col.name = "variant_id")
    if (!"variant_id" %in% names(tab)) stop("tab has no variant_id column", call. = FALSE)
    va <- variant_annotations(db, variant.ids = unique(tab$variant_id))
    va <- va[, c("variant_id", setdiff(names(va), names(tab))), drop = FALSE]
    keep_provenance(merge(tab, va, by = "variant_id", all.x = TRUE), tab)
}

#' @rdname add_context
#' @export
add_cnv_context <- function(tab, db = NULL) {
    tab <- ensure_long(tab, col.name = "cnv_id")
    if (!"cnv_id" %in% names(tab)) stop("tab has no cnv_id column", call. = FALSE)
    cn <- cnvs(db)
    cn <- cn[, c("cnv_id", setdiff(names(cn), names(tab))), drop = FALSE]
    keep_provenance(merge(tab, cn, by = "cnv_id", all.x = TRUE), tab)
}

#' Taxonomy helpers for microbiome measurements
#'
#' \code{deduplicate_taxonomy} makes taxon names unique at each level by
#' prefixing ancestors where the same name occurs under different parents;
#' \code{aggregate_taxa} sums measurements to each taxonomic level
#' (optionally normalizing each sample to proportions). Work for OTUs
#' (\code{otu_*} columns) and SGBs (\code{sgb_*} columns).
#'
#' @param taxa Output of \code{\link{otus}} or \code{\link{sgbs}}
#' @param na.value Name for missing levels
#' @param sep Separator for prefixed names
#' @param tab Long measurements with \code{sample_id}, an id column matching
#'   \code{taxa}, and \code{value}
#' @param id.col Taxon id column shared by \code{tab} and \code{taxa}
#' @param normalize Per-sample proportions
#' @param wide Return samples x taxa matrices
#' @return \code{deduplicate_taxonomy}: taxa table; \code{aggregate_taxa}:
#'   named list, one element per level
#' @name taxonomy
NULL

taxonomy_levels <- function(taxa) {
    lv <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    prefix <- if (any(startsWith(names(taxa), "otu_"))) "otu_" else "sgb_"
    intersect(paste0(prefix, lv), names(taxa))
}

#' @rdname taxonomy
#' @export
deduplicate_taxonomy <- function(taxa, na.value = "Unclassified", sep = "_") {
    levels <- taxonomy_levels(taxa)
    out <- taxa
    for (i in seq_along(levels)) {
        cur <- levels[i]
        out[[cur]][is.na(out[[cur]])] <- na.value
        taxa[[cur]][is.na(taxa[[cur]])] <- na.value
        uniq <- unique(taxa[, levels[1:i], drop = FALSE])
        dup <- uniq[duplicated(uniq[[cur]]), cur]
        w <- out[[cur]] %in% dup
        if (any(w))
            out[w, cur] <- apply(taxa[w, levels[1:i], drop = FALSE], 1, paste, collapse = sep)
    }
    out
}

#' @rdname taxonomy
#' @export
aggregate_taxa <- function(tab, taxa, id.col, normalize = TRUE, wide = FALSE, na.value = "Unclassified") {
    taxa <- deduplicate_taxonomy(taxa, na.value = na.value)
    levels <- taxonomy_levels(taxa)
    tab <- tab[!is.na(tab$value), , drop = FALSE]
    tab <- merge(tab, taxa, by = id.col)
    out <- list()
    for (cur in levels) {
        m <- stats::aggregate(tab$value, by = list(sample_id = tab$sample_id, taxon = tab[[cur]]), FUN = sum)
        names(m)[3] <- "value"
        if (normalize) m$value <- m$value / stats::ave(m$value, m$sample_id, FUN = sum)
        if (wide) {
            m <- to_matrix(m, col = "taxon", fun.aggregate = sum)
            m[is.na(m)] <- 0
        }
        out[[cur]] <- m
    }
    out
}
