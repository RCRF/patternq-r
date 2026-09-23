# Measurement-set level queries: which measurement sets exist, what they
# measure, and generic measurement retrieval.

#' Measurement sets with their assay and size
#'
#' @inheritParams samples
#' @return \code{assay_name}, \code{assay_technology},
#'   \code{measurement_set_name}, \code{measurement_count}
#' @export
measurement_sets <- function(db = NULL, ...) {
    df <- dataset_summary(db = db, ...)
    counts <- do_query(dq(find = list("?measurement-set-name", c("count", "?m")),
                          where = list(c("?ms", ":measurement-set/name", "?measurement-set-name"),
                                       c("?ms", ":measurement-set/measurements", "?m"))),
                       db = db, ...)
    df$measurement_count <- counts$count_m[match(df$measurement_set_name, counts$measurement_set_name)]
    df$measurement_count[is.na(df$measurement_count)] <- 0
    df
}

measurement_set_eid <- function(measurement.set, db) {
    eid <- do_query(dq(find = "?ms", `in` = list("?name"), args = list(measurement.set),
                       where = list(c("?ms", ":measurement-set/name", "?name"))),
                    db = db, simplify = TRUE)
    if (!length(eid)) stop(sprintf("No measurement set named '%s'", measurement.set), call. = FALSE)
    eid[1]
}

#' Measurement attributes in a measurement set
#'
#' \code{measurement_types} counts every attribute across the set's
#' measurements: value attributes (\code{tpm}, \code{percent-of-parent}, ...)
#' and target references (\code{gene-product}, \code{cell-population}, ...).
#' A set can mix several kinds of measurement, e.g. CyTOF populations with
#' \code{percent-of-parent} alongside population x marker
#' \code{median-channel-value}.
#'
#' \code{measurement_set_attributes} returns the attributes of a random
#' sample of \code{n} measurements, optionally only those carrying
#' \code{measurement}; cheap even for sets with millions of measurements.
#'
#' @inheritParams samples
#' @param measurement.set Measurement set name
#' @param measurement Optional measurement attribute (without namespace) the
#'   sampled measurements must carry
#' @param n Sample size
#' @return \code{measurement_types}: \code{attribute}, \code{kind}
#'   (\code{"value"} or \code{"target"}), \code{count};
#'   \code{measurement_set_attributes}: character vector of attribute idents
#' @export
measurement_types <- function(measurement.set, db = NULL, ...) {
    r <- do_query(dq(find = list("?attribute", c("count", "?m")),
                     `in` = list("?ms-name"), args = list(measurement.set),
                     where = list(c("?ms", ":measurement-set/name", "?ms-name"),
                                  c("?ms", ":measurement-set/measurements", "?m"),
                                  c("?m", "?a"), c("?a", ":db/ident", "?attribute"))),
                  db = db, timeout = 120, ...)
    r <- r[!r$attribute %in% c(":measurement/id", ":measurement/uid", ":measurement/sample"), , drop = FALSE]
    names(r)[names(r) == "count_m"] <- "count"
    r$kind <- ifelse(r$attribute %in% names(measurement_targets), "target", "value")
    r$attribute <- sub("^:measurement/", "", r$attribute)
    keep_provenance(order_columns(r[order(r$kind, -r$count), , drop = FALSE], c("attribute", "kind", "count")), r)
}

#' @rdname measurement_types
#' @export
measurement_set_attributes <- function(measurement.set, db = NULL, measurement = NULL, n = 200) {
    db <- ensure_db(db)
    where <- list(c("?ms", ":measurement-set/name", "?ms-name"),
                  c("?ms", ":measurement-set/measurements", "?m"))
    if (!is.null(measurement))
        where <- c(where, list(c("?m", paste0(":measurement/", sub("^:?measurement/", "", measurement)))))
    s <- raw_query(dq(find = list(list("sample", as.integer(n), "?m")), `in` = list("?ms-name"), args = list(measurement.set),
                      where = where), db = db, timeout = 120)
    eids <- unlist(s$query_result)
    if (!length(eids)) return(character(0))
    attrs <- do_query(dq(find = "?attr", `in` = list(c("?m", "...")), args = list(I(eids)),
                         where = list(c("?m", "?a"), c("?a", ":db/ident", "?attr"))),
                      db = db, simplify = TRUE)
    sort(setdiff(as.vector(attrs), c(":measurement/uid", ":measurement/id")))
}

# Measurement target references: the entity a measurement is "of", how to
# name it, and the result column.
measurement_targets <- list(
    ":measurement/gene-product" = list(clauses = list(c("?tgp", ":gene-product/gene", "?tg"),
                                                      c("?tg", ":gene/hgnc-symbol", "?hgnc-symbol")),
                                       ref = "?tgp", var = "?hgnc-symbol"),
    ":measurement/variant" = list(clauses = list(c("?tv", ":variant/id", "?variant-id")),
                                  ref = "?tv", var = "?variant-id"),
    ":measurement/cnv" = list(clauses = list(c("?tc", ":cnv/id", "?cnv-id")), ref = "?tc", var = "?cnv-id"),
    ":measurement/epitope" = list(clauses = list(c("?te", ":epitope/id", "?epitope-id")),
                                  ref = "?te", var = "?epitope-id"),
    ":measurement/cell-population" = list(clauses = list(c("?tcp", ":cell-population/name", "?cell-population")),
                                          ref = "?tcp", var = "?cell-population"),
    ":measurement/tcr" = list(clauses = list(c("?tt", ":tcr/id", "?tcr-id")), ref = "?tt", var = "?tcr-id"),
    ":measurement/otu" = list(clauses = list(c("?to", ":otu/id", "?otu-id")), ref = "?to", var = "?otu-id"),
    ":measurement/sgb" = list(clauses = list(c("?ts", ":sgb/metaphlan-id", "?sgb-id")), ref = "?ts", var = "?sgb-id"),
    ":measurement/pathway" = list(clauses = list(c("?tp", ":pathway/id", "?pathway-id")),
                                  ref = "?tp", var = "?pathway-id"),
    ":measurement/metabolite-feature" = list(clauses = list(c("?tmf", ":metabolite-feature/rt-mz-peak", "?metabolite-feature")),
                                             ref = "?tmf", var = "?metabolite-feature"),
    ":measurement/nanostring-signature" = list(clauses = list(c("?tns", ":nanostring-signature/name", "?signature")),
                                               ref = "?tns", var = "?signature"),
    ":measurement/atac-peak" = list(clauses = list(c("?tap", ":atac-peak/name", "?atac-peak")),
                                    ref = "?tap", var = "?atac-peak"),
    ":measurement/single-cell" = list(clauses = list(c("?tsc", ":single-cell/id", "?single-cell-id")),
                                      ref = "?tsc", var = "?single-cell-id")
)

enum_measurement_attrs <- c(":measurement/cnv-call", ":measurement/msi-status")

#' Measurement values of one type
#'
#' Generic retrieval of any measurement attribute (e.g. \code{"tpm"},
#' \code{"percent-of-parent"}, \code{"olink-npx"}, \code{"median-channel-value"})
#' with what each measurement is of (its target: gene, cell population,
#' epitope, variant, ...). Targets are detected from the measurement set's
#' measurements unless given.
#'
#' @inheritParams samples
#' @param measurement Measurement attribute without namespace
#' @param measurement.set Measurement set name. Required unless \code{targets}
#'   is given (targets are detected per measurement set).
#' @param samples Optional sample ids
#' @param targets Optional target reference attributes, e.g.
#'   \code{c(":measurement/cell-population", ":measurement/epitope")}; use
#'   \code{character(0)} for none
#' @param wide Return a samples x targets matrix (see \code{\link{to_matrix}})
#'   instead of long format
#' @param fun.aggregate Aggregation for duplicate sample/target cells when
#'   \code{wide = TRUE}
#' @return Long format: \code{sample_id}, \code{measurement_set}, target
#'   column(s), \code{value}; or a matrix when \code{wide = TRUE}
#' @export
measurements <- function(measurement, measurement.set = NULL, db = NULL, samples = NULL,
                         targets = NULL, wide = FALSE, fun.aggregate = mean, ...) {
    db <- ensure_db(db)
    if (is.null(targets)) {
        if (is.null(measurement.set))
            stop("Give measurement.set (or targets) so measurement targets can be detected", call. = FALSE)
        targets <- intersect(measurement_set_attributes(measurement.set, db = db, measurement = measurement),
                             names(measurement_targets))
    }
    q <- measurements_query(measurement, measurement.set = measurement.set, samples = samples,
                            targets = targets)
    df <- do_query(q, db = db, ...)
    if (is.character(df$value)) df$value <- ident_name(df$value)
    if (wide) {
        target.cols <- setdiff(names(df), c("sample_id", "measurement_set", "value"))
        if (length(target.cols) == 0) target.cols <- "measurement_set"
        return(keep_provenance(to_matrix(df, col = target.cols, fun.aggregate = fun.aggregate), df))
    }
    df
}

#' @rdname measurements
#' @export
measurements_query <- function(measurement, measurement.set = NULL, samples = NULL, targets = character(0)) {
    attr <- paste0(":measurement/", sub("^:?measurement/", "", measurement))
    value.var <- if (attr %in% enum_measurement_attrs) "?value-ref" else "?value"
    where <- list(c("?m", attr, value.var),
                  c("?m", ":measurement/sample", "?s"),
                  c("?s", ":sample/id", "?sample-id"),
                  c("?ms", ":measurement-set/measurements", "?m"),
                  c("?ms", ":measurement-set/name", "?measurement-set"))
    if (attr %in% enum_measurement_attrs) where <- c(where, list(c("?value-ref", ":db/ident", "?value")))
    find <- list("?sample-id", "?measurement-set")
    for (t in targets) {
        spec <- measurement_targets[[t]]
        if (is.null(spec)) stop("Unknown measurement target ", t, call. = FALSE)
        where <- c(where, list(c("?m", t, spec$ref)), spec$clauses)
        find <- c(find, list(spec$var))
    }
    find <- c(find, list("?value"))
    ins <- list(); args <- list()
    if (!is.null(measurement.set)) { ins <- c(ins, list("?measurement-set")); args <- c(args, list(measurement.set)) }
    if (!is.null(samples)) { ins <- c(ins, list(c("?sample-id", "..."))); args <- c(args, list(I(samples))) }
    q <- dq(find = find, with = "?m", where = where, `in` = if (length(ins)) ins else NULL, args = args)
    q
}

#' Which samples were measured in which measurement sets
#'
#' @inheritParams samples
#' @return \code{subject_id}, \code{sample_id}, \code{timepoint_id},
#'   \code{assay_name}, \code{assay_technology}, \code{measurement_set_name}
#' @export
sample_assays <- function(db = NULL, ...) {
    # measurement set x sample, per measurement set, so the query only walks
    # each set's measurement -> sample refs
    db <- ensure_db(db)
    sets <- dataset_summary(db)
    parts <- lapply(seq_len(nrow(sets)), function(i) {
        ms <- sets$measurement_set_name[i]
        r <- do_query(dq(find = "?sample-id", `in` = list("?ms-name"), args = list(ms),
                         where = list(c("?ms", ":measurement-set/name", "?ms-name"),
                                      c("?ms", ":measurement-set/measurements", "?m"),
                                      c("?m", ":measurement/sample", "?s"),
                                      c("?s", ":sample/id", "?sample-id"))),
                      db = db, simplify = TRUE, ...)
        if (!length(r)) return(NULL)
        data.frame(sample_id = r, assay_name = sets$assay_name[i],
                   assay_technology = sets$assay_technology[i],
                   measurement_set_name = ms, stringsAsFactors = FALSE)
    })
    df <- bind_rows_fill(parts)
    if (!nrow(df)) return(df)
    smp <- samples(db)
    keep <- intersect(c("sample_id", "subject_id", "timepoint_id"), names(smp))
    df <- merge(smp[, keep, drop = FALSE], df, by = "sample_id", all.y = TRUE)
    order_columns(df, c("subject_id", "sample_id", "timepoint_id"))
}

#' Measurement matrices (file-backed measurements)
#'
#' @inheritParams samples
#' @return \code{assay_name}, \code{measurement_set_name},
#'   \code{matrix_name}, \code{measurement_type}, \code{matrix_key}; pass
#'   \code{matrix_key} to \code{\link{measurement_matrix}}
#' @export
measurement_matrices <- function(db = NULL, ...) {
    df <- do_query(measurement_matrices_query(), db = db, ...)
    df$measurement_type <- ident_name(df$measurement_type)
    df
}

#' @rdname measurement_matrices
#' @export
measurement_matrices_query <- function() {
    dq(find = c("?assay-name", "?measurement-set-name", "?matrix-name", "?measurement-type", "?matrix-key"),
       where = list(c("?a", ":assay/name", "?assay-name"),
                    c("?a", ":assay/measurement-sets", "?ms"),
                    c("?ms", ":measurement-set/name", "?measurement-set-name"),
                    c("?ms", ":measurement-set/measurement-matrices", "?mm"),
                    c("?mm", ":measurement-matrix/name", "?matrix-name"),
                    c("?mm", ":measurement-matrix/measurement-type", "?mt"),
                    c("?mt", ":db/ident", "?measurement-type"),
                    c("?mm", ":measurement-matrix/backing-file", "?matrix-key")))
}

#' Download a measurement matrix by name
#'
#' @inheritParams samples
#' @param matrix.name Measurement matrix name (see \code{\link{measurement_matrices}})
#' @return \code{data.frame}
#' @export
measurement_matrix_by_name <- function(matrix.name, db = NULL) {
    mm <- measurement_matrices(db)
    key <- mm$matrix_key[mm$matrix_name == matrix.name]
    if (!length(key)) stop(sprintf("No measurement matrix named '%s'", matrix.name), call. = FALSE)
    measurement_matrix(key[1], db = db)
}

#' Isoform-level expression for a gene
#'
#' @inheritParams samples
#' @param gene HGNC symbol
#' @param samples Optional sample ids
#' @return \code{sample_id}, \code{transcript_id}, \code{transcript_length},
#'   \code{isoform_percent}, \code{effective_length}
#' @export
isoforms <- function(gene, db = NULL, samples = NULL, ...) {
    q <- dq(find = c("?sample-id", "?transcript-id", "?transcript-length", "?isoform-percent", "?effective-length"),
            with = "?m",
            `in` = c(list("?hgnc"), if (!is.null(samples)) list(c("?sample-id", "..."))),
            args = c(list(gene), if (!is.null(samples)) list(I(samples))),
            where = list(c("?g", ":gene/hgnc-symbol", "?hgnc"),
                         c("?gp", ":gene-product/gene", "?g"),
                         c("?gp", ":gene-product/id", "?transcript-id"),
                         c("?gp", ":gene-product/transcript-length", "?transcript-length"),
                         c("?m", ":measurement/gene-product", "?gp"),
                         c("?m", ":measurement/isoform-percent", "?isoform-percent"),
                         c("?m", ":measurement/effective-transcript-length", "?effective-length"),
                         c("?m", ":measurement/sample", "?s"),
                         c("?s", ":sample/id", "?sample-id")))
    do_query(q, db = db, ...)
}

# CNV data is large (segments x samples, or genes x samples); the CNV queries
# insist on a subset. Same rule in the Python and Clojure libraries.
require_cnv_subset <- function(genes, samples, subjects) {
    if (is.null(genes) && is.null(samples) && is.null(subjects))
        stop("CNV queries need a subset: give genes, samples or subjects", call. = FALSE)
}

subject_clauses <- function(subjects) {
    if (is.null(subjects)) return(list(where = list(), ins = list(), args = list()))
    list(where = list(c("?s", ":sample/subject", "?p"), c("?p", ":subject/id", "?subject-id")),
         ins = list(c("?subject-id", "...")), args = list(I(subjects)))
}

#' Copy number segments per sample
#'
#' Segment-level CNV measurements (e.g. CNVkit / ASCAT calls) with the genes
#' each segment overlaps.
#'
#' @inheritParams samples
#' @param genes HGNC symbols; restricts to segments overlapping them and
#'   returns one row per segment and gene
#' @param samples Sample ids
#' @param subjects Subject ids. At least one of \code{genes}, \code{samples},
#'   \code{subjects} is required.
#' @return \code{sample_id}, \code{measurement_set}, \code{cnv_id},
#'   \code{contig}, \code{start}, \code{end}, \code{segment_mean_lrr},
#'   \code{absolute_cn} (where present), \code{hgnc_symbol} (with \code{genes})
#' @export
cnv_segments <- function(db = NULL, genes = NULL, samples = NULL, subjects = NULL, ...) {
    require_cnv_subset(genes, samples, subjects)
    where <- list(c("?m", ":measurement/cnv", "?c"),
                  c("?c", ":cnv/id", "?cnv-id"),
                  c("?m", ":measurement/sample", "?s"),
                  c("?s", ":sample/id", "?sample-id"),
                  c("?ms", ":measurement-set/measurements", "?m"),
                  c("?ms", ":measurement-set/name", "?measurement-set"))
    find <- list("?sample-id", "?measurement-set", "?cnv-id",
                 list("pull", "?c", list(list(":cnv/genomic-coordinates" = c(":genomic-coordinate/contig",
                                                                          ":genomic-coordinate/start",
                                                                          ":genomic-coordinate/end")))),
                 list("pull", "?m", c(":measurement/segment-mean-lrr", ":measurement/absolute-cn",
                                      ":measurement/a-allele-cn", ":measurement/b-allele-cn",
                                      ":measurement/loh")))
    ins <- list(); args <- list()
    if (!is.null(genes)) {
        where <- c(where, list(c("?c", ":cnv/genes", "?g"), c("?g", ":gene/hgnc-symbol", "?hgnc-symbol")))
        find <- c(find, list("?hgnc-symbol"))
        ins <- c(ins, list(c("?hgnc-symbol", "..."))); args <- c(args, list(I(genes)))
    }
    if (!is.null(samples)) { ins <- c(ins, list(c("?sample-id", "..."))); args <- c(args, list(I(samples))) }
    sc <- subject_clauses(subjects)
    where <- c(where, sc$where); ins <- c(ins, sc$ins); args <- c(args, sc$args)
    df <- do_query(dq(find = find, where = where, `in` = if (length(ins)) ins else NULL, args = args),
                   db = db, ...)
    names(df) <- sub("^(genomic_coordinate|measurement)_(?!set$)", "", names(df), perl = TRUE)
    keep_provenance(order_columns(df, c("sample_id", "measurement_set", "cnv_id", "hgnc_symbol",
                                        "contig", "start", "end")), df)
}

#' Gene-level copy number calls
#'
#' Gene-level CNV measurements (e.g. GISTIC2 discrete calls,
#' \code{:measurement/cnv-call-score}, or \code{:measurement/cnv-call}).
#'
#' @inheritParams samples
#' @param genes HGNC symbols
#' @param samples Sample ids
#' @param subjects Subject ids. At least one of \code{genes}, \code{samples},
#'   \code{subjects} is required.
#' @param measurement \code{"cnv-call-score"} (default) or \code{"cnv-call"}
#' @return \code{sample_id}, \code{measurement_set}, \code{hgnc_symbol}, \code{value}
#' @export
cnv_gene_calls <- function(db = NULL, genes = NULL, samples = NULL, subjects = NULL,
                           measurement = "cnv-call-score", ...) {
    require_cnv_subset(genes, samples, subjects)
    attr <- paste0(":measurement/", sub("^:?measurement/", "", measurement))
    enum <- attr %in% enum_measurement_attrs
    # gene-first clause order: gene-level CNV sets are large (genes x samples)
    where <- list(c("?g", ":gene/hgnc-symbol", "?hgnc-symbol"),
                  c("?gp", ":gene-product/gene", "?g"),
                  c("?m", ":measurement/gene-product", "?gp"),
                  c("?m", attr, if (enum) "?value-ref" else "?value"),
                  c("?m", ":measurement/sample", "?s"),
                  c("?s", ":sample/id", "?sample-id"),
                  c("?ms", ":measurement-set/measurements", "?m"),
                  c("?ms", ":measurement-set/name", "?measurement-set"))
    if (enum) where <- c(where, list(c("?value-ref", ":db/ident", "?value")))
    ins <- list(); args <- list()
    if (!is.null(genes)) { ins <- c(ins, list(c("?hgnc-symbol", "..."))); args <- c(args, list(I(genes))) }
    if (!is.null(samples)) { ins <- c(ins, list(c("?sample-id", "..."))); args <- c(args, list(I(samples))) }
    sc <- subject_clauses(subjects)
    where <- c(where, sc$where); ins <- c(ins, sc$ins); args <- c(args, sc$args)
    df <- do_query(dq(find = c("?sample-id", "?measurement-set", "?hgnc-symbol", "?value"), with = "?m",
                      where = where, `in` = if (length(ins)) ins else NULL, args = args),
                   db = db, ...)
    if (is.character(df$value)) df$value <- ident_name(df$value)
    df
}
