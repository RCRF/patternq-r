# Canned queries over a dataset database.
#
# Every dataset is its own database, so queries start from samples,
# subjects, assays and measurement sets directly; there is no dataset
# argument. Each function has a *_query() companion returning the query as
# data, so it can be inspected, modified (c_query) or run elsewhere.

#' Samples
#'
#' @param db Database name (defaults to \code{\link{current_db}})
#' @param ... Passed to \code{\link{do_query}}
#' @return One row per sample: \code{sample_id}, \code{subject_id},
#'   \code{timepoint_id}, and the other sample attributes present
#' @export
samples <- function(db = NULL, ...) {
    do_query(samples_query(), db = db, ...)
}

#' @rdname samples
#' @export
samples_query <- function() {
    dq(find = list(list("pull", "?s", list(
        "*",
        list(":sample/subject" = ":subject/id"),
        list(":sample/timepoint" = c(":timepoint/id", ":timepoint/relative-order")),
        list(":sample/specimen" = ":db/ident"),
        list(":sample/type" = ":db/ident"),
        list(":sample/container" = ":db/ident"),
        list(":sample/study-day" = ":study-day/id"),
        list(":sample/gdc-anatomic-site" = ":gdc-anatomic-site/name")))),
       where = list(c("?s", ":sample/id")))
}

#' Subjects (participants)
#'
#' @inheritParams samples
#' @return One row per subject: \code{subject_id} and demographic attributes
#'   present (enums as names, e.g. \code{subject_sex = "female"})
#' @export
subjects <- function(db = NULL, ...) {
    df <- do_query(subjects_query(), db = db, ...)
    df
}

#' @rdname subjects
#' @export
subjects_query <- function() {
    dq(find = list(list("pull", "?s", list(
        "*",
        list(":subject/sex" = ":db/ident"),
        list(":subject/race" = ":db/ident"),
        list(":subject/ethnicity" = ":db/ident"),
        list(":subject/smoker" = ":db/ident"),
        list(":subject/cause-of-death" = ":db/ident"),
        list(":subject/meddra-disease" = ":meddra-disease/preferred-name")))),
       where = list(c("?s", ":subject/id")))
}

#' Assays and measurement sets in a dataset
#'
#' @inheritParams samples
#' @return One row per measurement set: \code{assay_name},
#'   \code{assay_technology}, \code{measurement_set_name}
#' @export
dataset_summary <- function(db = NULL, ...) {
    df <- do_query(dataset_summary_query(), db = db, ...)
    df$assay_technology <- ident_name(df$assay_technology)
    out <- df[order(df$assay_name, df$measurement_set_name), , drop = FALSE]
    rownames(out) <- NULL
    keep_provenance(out, df)
}

#' @rdname dataset_summary
#' @export
dataset_summary_query <- function() {
    dq(find = c("?assay-name", "?assay-technology", "?measurement-set-name"),
       where = list(c("?a", ":assay/name", "?assay-name"),
                    c("?a", ":assay/technology", "?t"),
                    c("?t", ":db/ident", "?assay-technology"),
                    c("?a", ":assay/measurement-sets", "?ms"),
                    c("?ms", ":measurement-set/name", "?measurement-set-name")))
}

#' Somatic variant measurements
#'
#' @inheritParams samples
#' @param samples Optional character vector of sample ids
#' @param genes Optional character vector of HGNC symbols
#' @param measurement.set Optional measurement set name (datasets can have
#'   several variant measurement sets, e.g. WES and WGS)
#' @return One row per variant measurement: \code{sample_id},
#'   \code{measurement_set}, \code{variant_id}, \code{hgnc_symbol},
#'   \code{vaf}, plus \code{HGVSp}, \code{impact}, \code{t_depth} where
#'   present
#' @export
variants <- function(db = NULL, samples = NULL, genes = NULL,
                     measurement.set = NULL, ...) {
    df <- do_query(variants_query(samples = samples, genes = genes,
                                  measurement.set = measurement.set),
                   db = db, ...)
    names(df)[names(df) == "variant_HGVSp"] <- "HGVSp"
    names(df)[names(df) == "variant_HGVSc"] <- "HGVSc"
    names(df)[names(df) == "variant_impact"] <- "impact"
    names(df)[names(df) == "gene_hgnc_symbol"] <- "hgnc_symbol"
    names(df)[names(df) == "measurement_t_depth"] <- "t_depth"
    # HGVSp / HGVSc are cardinality many; collapse to "; "-joined strings
    for (col in intersect(c("HGVSp", "HGVSc"), names(df))) df[[col]] <- join_many(df[[col]])
    order_columns(df, c("sample_id", "measurement_set", "variant_id", "hgnc_symbol",
                        "HGVSp", "HGVSc", "impact", "vaf", "t_depth"))
}

#' @rdname variants
#' @export
variants_query <- function(samples = NULL, genes = NULL, measurement.set = NULL) {
    where <- list(c("?m", ":measurement/vaf", "?vaf"),
                  c("?m", ":measurement/variant", "?v"),
                  c("?m", ":measurement/sample", "?s"),
                  c("?s", ":sample/id", "?sample-id"),
                  c("?ms", ":measurement-set/measurements", "?m"),
                  c("?ms", ":measurement-set/name", "?measurement-set"))
    ins <- list()
    args <- list()
    if (!is.null(samples)) {
        ins <- c(ins, list(c("?sample-id", "...")))
        args <- c(args, list(I(samples)))
    }
    if (!is.null(genes)) {
        where <- c(where, list(c("?v", ":variant/gene", "?g"),
                               c("?g", ":gene/hgnc-symbol", "?gene")))
        ins <- c(ins, list(c("?gene", "...")))
        args <- c(args, list(I(genes)))
    }
    if (!is.null(measurement.set)) {
        ins <- c(ins, list("?measurement-set"))
        args <- c(args, list(measurement.set))
    }
    dq(find = list("?sample-id", "?measurement-set", "?vaf",
                   list("pull", "?v", list(":variant/id", ":variant/HGVSp", ":variant/HGVSc",
                                           list(":variant/gene" = ":gene/hgnc-symbol"),
                                           list(":variant/impact" = ":db/ident"))),
                   list("pull", "?m", list(":measurement/t-depth"))),
       where = where, `in` = if (length(ins)) ins else NULL, args = args)
}

rnaseq_attrs <- c("tpm", "fpkm", "fpkm-upper-quartile", "rsem-normalized-count",
                  "rsem-raw-count", "rsem-scaled-estimate", "read-count", "rpkm",
                  "kallisto-abundance")

#' Gene expression measurements
#'
#' @inheritParams samples
#' @param genes Character vector of HGNC symbols (\code{NULL} for all genes;
#'   can be large)
#' @param samples Optional character vector of sample ids
#' @param measurement Measurement attribute without namespace, e.g.
#'   \code{"tpm"}, \code{"rsem-normalized-count"}, \code{"read-count"}
#' @param measurement.set Optional measurement set name. Datasets may have
#'   several (bulk RNA-seq, pseudobulk, ...); without it values from every set
#'   carrying the attribute are returned and \code{measurement_set} tells them
#'   apart.
#' @return Long format: \code{sample_id}, \code{hgnc_symbol},
#'   \code{measurement_set}, \code{value}
#' @export
gene_expression <- function(db = NULL, genes = NULL, samples = NULL,
                            measurement = "tpm", measurement.set = NULL, ...) {
    df <- do_query(gene_expression_query(genes = genes, samples = samples,
                                         measurement = measurement,
                                         measurement.set = measurement.set),
                   db = db, ...)
    df
}

#' @rdname gene_expression
#' @export
gene_expression_query <- function(genes = NULL, samples = NULL, measurement = "tpm",
                                  measurement.set = NULL) {
    attr <- paste0(":measurement/", sub("^:?measurement/", "", measurement))
    where <- list(c("?g", ":gene/hgnc-symbol", "?hgnc-symbol"),
                  c("?gp", ":gene-product/gene", "?g"),
                  c("?m", ":measurement/gene-product", "?gp"),
                  c("?m", attr, "?value"),
                  c("?m", ":measurement/sample", "?s"),
                  c("?s", ":sample/id", "?sample-id"),
                  c("?ms", ":measurement-set/measurements", "?m"),
                  c("?ms", ":measurement-set/name", "?measurement-set"))
    ins <- list()
    args <- list()
    if (!is.null(genes)) {
        ins <- c(ins, list(c("?hgnc-symbol", "...")))
        args <- c(args, list(I(genes)))
    }
    if (!is.null(samples)) {
        ins <- c(ins, list(c("?sample-id", "...")))
        args <- c(args, list(I(samples)))
    }
    if (!is.null(measurement.set)) {
        ins <- c(ins, list("?measurement-set"))
        args <- c(args, list(measurement.set))
    }
    dq(find = c("?sample-id", "?hgnc-symbol", "?measurement-set", "?value"),
       with = "?m", where = where, `in` = if (length(ins)) ins else NULL, args = args)
}
