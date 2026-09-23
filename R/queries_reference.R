# Reference data queries. Reference entities (genes, gene products, proteins,
# epitopes, cell types, ...) are included in each dataset database, so these
# take a db like every other query.

run_simple <- function(q, db, simplify = FALSE, ...) {
    do_query(q, db = db, simplify = simplify, ...)
}

#' HGNC gene symbols
#' @inheritParams samples
#' @return Character vector
#' @export
gene_symbols <- function(db = NULL, ...) {
    run_simple(gene_symbols_query(), db, simplify = TRUE, ...)
}

#' @rdname gene_symbols
#' @export
gene_symbols_query <- function() {
    dq(find = "?hgnc-symbol", where = list(c("_", ":gene/hgnc-symbol", "?hgnc-symbol")))
}

#' Genes
#'
#' @inheritParams samples
#' @return One row per gene: \code{gene_hgnc_symbol}, \code{gene_hgnc_name},
#'   ids, and list columns \code{gene_previous_hgnc_symbols},
#'   \code{gene_alias_hgnc_symbols}
#' @export
genes <- function(db = NULL, ...) {
    run_simple(genes_query(), db, ...)
}

#' @rdname genes
#' @export
genes_query <- function() {
    dq(find = list(list("pull", "?g", list(":gene/hgnc-symbol", ":gene/hgnc-id", ":gene/hgnc-name",
                                           ":gene/ensembl-id", ":gene/previous-hgnc-symbols",
                                           ":gene/alias-hgnc-symbols",
                                           list(":gene/hgnc-locus-group" = ":db/ident")))),
       where = list(c("?g", ":gene/hgnc-symbol")))
}

#' Gene products and their genes
#' @inheritParams samples
#' @return \code{gene_product_id}, \code{hgnc_symbol}
#' @export
gene_products <- function(db = NULL, ...) {
    run_simple(gene_products_query(), db, ...)
}

#' @rdname gene_products
#' @export
gene_products_query <- function() {
    dq(find = c("?gene-product-id", "?hgnc-symbol"),
       where = list(c("?gp", ":gene-product/id", "?gene-product-id"),
                    c("?gp", ":gene-product/gene", "?g"),
                    c("?g", ":gene/hgnc-symbol", "?hgnc-symbol")))
}

#' Genes with genomic coordinates
#' @inheritParams samples
#' @param genes Optional HGNC symbols
#' @return \code{hgnc_symbol}, \code{assembly}, \code{contig}, \code{strand},
#'   \code{start}, \code{end}
#' @export
gene_coordinates <- function(db = NULL, genes = NULL, ...) {
    df <- run_simple(gene_coordinates_query(genes), db, ...)
    df$assembly <- ident_name(df$assembly)
    df
}

#' @rdname gene_coordinates
#' @export
gene_coordinates_query <- function(genes = NULL) {
    q <- dq(find = c("?hgnc-symbol", "?assembly", "?contig", "?strand", "?start", "?end"),
            where = list(c("?g", ":gene/hgnc-symbol", "?hgnc-symbol"),
                         c("?g", ":gene/genomic-coordinates", "?gc"),
                         c("?gc", ":genomic-coordinate/assembly", "?a"),
                         c("?a", ":db/ident", "?assembly"),
                         c("?gc", ":genomic-coordinate/contig", "?contig"),
                         c("?gc", ":genomic-coordinate/strand", "?strand"),
                         c("?gc", ":genomic-coordinate/start", "?start"),
                         c("?gc", ":genomic-coordinate/end", "?end")))
    if (!is.null(genes)) {
        q$query[["in"]] <- list(c("?hgnc-symbol", "..."))
        q$args <- list(I(genes))
    }
    q
}

#' Variant annotations
#'
#' Variant reference entities (not measurements): id, gene, HGVS, impact,
#' consequences, classification.
#'
#' @inheritParams samples
#' @param variant.ids Optional variant ids
#' @param genes Optional HGNC symbols
#' @return One row per variant
#' @export
variant_annotations <- function(db = NULL, variant.ids = NULL, genes = NULL, ...) {
    df <- run_simple(variant_annotations_query(variant.ids, genes), db, ...)
    for (col in intersect(c("variant_HGVSp", "variant_HGVSc", "variant_so_consequences", "variant_external_ids"), names(df)))
        df[[col]] <- join_many(df[[col]])
    names(df) <- sub("^variant_", "", names(df))
    names(df)[names(df) == "gene_hgnc_symbol"] <- "hgnc_symbol"
    names(df)[names(df) == "id"] <- "variant_id"
    keep_provenance(order_columns(df, c("variant_id", "hgnc_symbol", "HGVSp", "HGVSc", "impact")), df)
}

#' @rdname variant_annotations
#' @export
variant_annotations_query <- function(variant.ids = NULL, genes = NULL) {
    where <- list(c("?v", ":variant/id", "?variant-id"))
    ins <- list(); args <- list()
    if (!is.null(variant.ids)) {
        ins <- c(ins, list(c("?variant-id", "..."))); args <- c(args, list(I(variant.ids)))
    }
    if (!is.null(genes)) {
        where <- c(where, list(c("?v", ":variant/gene", "?g"), c("?g", ":gene/hgnc-symbol", "?gene")))
        ins <- c(ins, list(c("?gene", "..."))); args <- c(args, list(I(genes)))
    }
    dq(find = list(list("pull", "?v", list(
        ":variant/id", ":variant/HGVSp", ":variant/HGVSc", ":variant/ref-allele", ":variant/alt-allele",
        ":variant/coordinate-string", ":variant/dbSNP", ":variant/max-af", ":variant/external-ids",
        list(":variant/gene" = ":gene/hgnc-symbol"),
        list(":variant/impact" = ":db/ident"),
        list(":variant/classification" = ":db/ident"),
        list(":variant/type" = ":db/ident"),
        list(":variant/so-consequences" = ":so-sequence-feature/name")))),
       where = where, `in` = if (length(ins)) ins else NULL, args = args)
}

#' CNV segments (reference entities)
#' @inheritParams samples
#' @return One row per CNV: \code{cnv_id}, coordinates, and list column
#'   \code{genes}
#' @export
cnvs <- function(db = NULL, ...) {
    run_simple(cnvs_query(), db, ...)
}

#' @rdname cnvs
#' @export
cnvs_query <- function() {
    dq(find = list(list("pull", "?c", list(
        ":cnv/id",
        list(":cnv/genomic-coordinates" = c(":genomic-coordinate/contig", ":genomic-coordinate/start",
                                             ":genomic-coordinate/end")),
        list(":cnv/genes" = ":gene/hgnc-symbol")))),
       where = list(c("?c", ":cnv/id")))
}

simple_name_query <- function(attr, var) {
    dq(find = var, where = list(c("_", attr, var)))
}

#' Reference vocabularies
#'
#' Names of reference entities present in a database.
#'
#' @inheritParams samples
#' @return Character vector
#' @name reference_names
NULL

#' @rdname reference_names
#' @export
gdc_anatomic_sites <- function(db = NULL, ...)
    run_simple(simple_name_query(":gdc-anatomic-site/name", "?name"), db, simplify = TRUE, ...)

#' @rdname reference_names
#' @export
proteins <- function(db = NULL, ...)
    run_simple(dq(find = list(list("pull", "?p", list(":protein/preferred-name", ":protein/uniprot-name",
                                                      ":protein/uniprot-accessions",
                                                      list(":protein/gene" = ":gene/hgnc-symbol")))),
                  where = list(c("?p", ":protein/preferred-name"))), db, ...)

#' @rdname reference_names
#' @export
epitopes <- function(db = NULL, ...)
    run_simple(simple_name_query(":epitope/id", "?epitope-id"), db, simplify = TRUE, ...)

#' @rdname reference_names
#' @export
cell_types <- function(db = NULL, ...)
    run_simple(simple_name_query(":cell-type/co-name", "?co-name"), db, simplify = TRUE, ...)

#' @rdname reference_names
#' @export
meddra_diseases <- function(db = NULL, ...)
    run_simple(simple_name_query(":meddra-disease/preferred-name", "?name"), db, simplify = TRUE, ...)

#' @rdname reference_names
#' @export
drugs <- function(db = NULL, ...)
    run_simple(simple_name_query(":drug/preferred-name", "?name"), db, simplify = TRUE, ...)

#' Dataset entity metadata
#'
#' The dataset entity stored in the database: name, description, doi, url.
#' @inheritParams samples
#' @return One-row data frame
#' @export
dataset_info <- function(db = NULL, ...) {
    run_simple(dq(find = list(list("pull", "?d", c(":dataset/name", ":dataset/description",
                                                   ":dataset/doi", ":dataset/url"))),
                  where = list(c("?d", ":dataset/name"))), db, ...)
}

#' Schema name and version of a database
#' @inheritParams samples
#' @return List with \code{name} and \code{version}
#' @export
schema_info <- function(db = NULL, ...) {
    r <- do_query(dq(find = c("?name", "?version"),
                     where = list(c("?e", ":unify.schema/version", "?version"),
                                  c("?e", ":unify.schema/name", "?name"))), db = db, ...)
    list(name = r$name[1], version = r$version[1])
}

#' All idents in a database
#' @inheritParams samples
#' @return \code{db_id}, \code{db_ident}
#' @export
db_idents <- function(db = NULL, ...) {
    run_simple(dq(find = c("?db-id", "?db-ident"), where = list(c("?db-id", ":db/ident", "?db-ident"))), db, ...)
}

#' Map gene symbols to current HGNC symbols
#'
#' Resolves previous and alias symbols to current HGNC symbols, case
#' insensitively.
#'
#' @param symbols Character vector of gene symbols
#' @param all.genes Output of \code{\link{genes}} (fetched from \code{db} if
#'   \code{NULL})
#' @param db Database name, used if \code{all.genes} is \code{NULL}
#' @param warn Warn about symbols that could not be mapped
#' @return Named character vector: names are the input symbols, values the
#'   current HGNC symbols (\code{NA} if unmapped)
#' @export
map_gene_symbols <- function(symbols, all.genes = NULL, db = NULL, warn = TRUE) {
    if (is.null(all.genes)) all.genes <- genes(db)
    lookup <- stats::setNames(all.genes$gene_hgnc_symbol, toupper(all.genes$gene_hgnc_symbol))
    for (col in c("gene_previous_hgnc_symbols", "gene_alias_hgnc_symbols")) {
        if (!col %in% names(all.genes)) next
        vals <- all.genes[[col]]
        for (i in seq_along(vals)) {
            v <- vals[[i]]
            if (length(v) == 0 || all(is.na(v))) next
            new <- toupper(unlist(v))
            new <- new[!new %in% names(lookup)]  # current symbols take precedence
            lookup <- c(lookup, stats::setNames(rep(all.genes$gene_hgnc_symbol[i], length(new)), new))
        }
    }
    out <- unname(lookup[toupper(symbols)])
    names(out) <- symbols
    if (warn && any(is.na(out)))
        warning(sprintf("%d of %d gene symbols could not be mapped", sum(is.na(out)), length(out)))
    out
}
