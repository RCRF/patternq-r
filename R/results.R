# Converting query results to data frames.

json_null_to_na <- function(v) if (is.null(v)) NA else v

column_from_rows <- function(rows, j) {
    vals <- lapply(rows, function(r) json_null_to_na(r[[j]]))
    scalar <- vapply(vals, function(v) is.atomic(v) && length(v) == 1, logical(1))
    if (all(scalar)) return(unlist(vals, use.names = FALSE))
    I(vals)
}

relation_to_df <- function(rows, col.names) {
    if (!length(rows)) {
        df <- as.data.frame(stats::setNames(replicate(length(col.names), logical(0), simplify = FALSE),
                                     col.names))
        return(df)
    }
    cols <- lapply(seq_along(col.names), function(j) column_from_rows(rows, j))
    names(cols) <- col.names
    df <- as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- col.names
    df
}

#' Keyword name without namespace
#'
#' \code{":variant.impact/high"} becomes \code{"high"}. Non-keyword values are
#' returned unchanged.
#'
#' @param x Character vector
#' @return Character vector
#' @export
ident_name <- function(x) {
    if (!is.character(x)) return(x)
    ifelse(!is.na(x) & grepl("^:[^/]+/", x), sub("^:[^/]+/", "", x), x)
}

is_map <- function(x) is.list(x) && !is.null(names(x)) && length(x) > 0

is_ident_map <- function(x) is_map(x) && identical(names(x), ":db/ident")

attr_col <- function(k) clean_names(k)

#' Flatten one pulled entity into a named list of columns
#'
#' Scalar attributes become columns named after the attribute
#' (\code{:subject/id} -> \code{subject_id}). Enum references
#' (\code{{:db/ident ...}}) become the enum name without namespace, in a column
#' named after the referring attribute. Nested single references are flattened
#' recursively. Cardinality-many values become list columns.
#'
#' @param m A pulled entity (named list from JSON)
#' @param exclude.ids Drop \code{db/id} and \code{*/uid} attributes
#' @return Named list of column values
#' @export
flatten_pull <- function(m, exclude.ids = TRUE) {
    out <- list()
    for (k in names(m)) {
        v <- m[[k]]
        if (exclude.ids && (k == ":db/id" || grepl("/uid$", k))) next
        if (is.null(v)) next
        if (is_ident_map(v)) {
            out[[attr_col(k)]] <- ident_name(v[[":db/ident"]])
        } else if (is_eid_map(v)) {
            if (!exclude.ids) out[[attr_col(k)]] <- v[[":db/id"]]
        } else if (is_map(v)) {
            nested <- flatten_pull(v, exclude.ids = exclude.ids)
            for (nk in names(nested)) {
                col <- if (is.null(out[[nk]])) nk else paste(attr_col(k), nk, sep = "_")
                out[[col]] <- nested[[nk]]
            }
        } else if (is.list(v)) {
            if (all(vapply(v, is_ident_map, logical(1)))) {
                out[[attr_col(k)]] <- list(vapply(v, function(e) ident_name(e[[":db/ident"]]), character(1)))
            } else if (all(vapply(v, function(e) is.atomic(e) && length(e) == 1, logical(1)))) {
                out[[attr_col(k)]] <- list(unlist(v))
            } else {
                out[[attr_col(k)]] <- list(bind_rows_fill(lapply(v, function(e) {
                    if (is_map(e)) as.data.frame(lapply(flatten_pull(e, exclude.ids), scalar_or_list),
                                                 stringsAsFactors = FALSE, check.names = FALSE)
                    else data.frame(value = e)
                })))
            }
        } else {
            out[[attr_col(k)]] <- v
        }
    }
    out
}

scalar_or_list <- function(v) if (is.list(v)) I(v) else v

pull_rows_to_df <- function(rows, q, exclude.ids = TRUE) {
    find <- q$query$find
    is_pull <- vapply(find, function(e) identical(as.character(unlist(e)[1]), "pull"), logical(1))
    other.names <- find_column_names(q)
    recs <- lapply(rows, function(r) {
        rec <- list()
        for (j in seq_along(find)) {
            v <- r[[j]]
            if (is_pull[j]) {
                flat <- flatten_pull(v, exclude.ids = exclude.ids)
                for (nk in names(flat)) rec[[nk]] <- flat[[nk]]
            } else {
                rec[[other.names[j]]] <- json_null_to_na(v)
            }
        }
        rec
    })
    records_to_df(recs)
}

# Combine a list of named lists (records) into a data frame, filling
# missing fields with NA; fields holding lists become list columns.
records_to_df <- function(recs) {
    if (!length(recs)) return(data.frame())
    cols <- unique(unlist(lapply(recs, names)))
    out <- lapply(cols, function(k) {
        vals <- lapply(recs, function(r) {
            v <- r[[k]]
            if (is.null(v)) NA else v
        })
        scalar <- vapply(vals, function(v) is.atomic(v) && length(v) == 1, logical(1))
        if (all(scalar)) unlist(vals, use.names = FALSE)
        else I(lapply(vals, function(v) if (is.list(v) && length(v) == 1 && !is.data.frame(v)) v[[1]] else v))
    })
    names(out) <- cols
    df <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- cols
    df
}

# rbind data frames with differing columns
bind_rows_fill <- function(dfs) {
    dfs <- dfs[!vapply(dfs, is.null, logical(1))]
    if (!length(dfs)) return(data.frame())
    cols <- unique(unlist(lapply(dfs, names)))
    dfs <- lapply(dfs, function(d) {
        for (k in setdiff(cols, names(d))) d[[k]] <- rep(NA, nrow(d))
        d[, cols, drop = FALSE]
    })
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    out
}

#' Provenance of a query result
#'
#' Results returned by patternq carry the database name, the database basis t
#' the query ran against, and the client time of the query, so an analysis can
#' be reproduced against the same database state.
#'
#' @param x A result returned by a patternq query function
#' @return A list with \code{db}, \code{basis_t} and \code{timestamp}, or
#'   \code{NULL}
#' @export
provenance <- function(x) {
    attr(x, "patternq_provenance")
}

with_provenance <- function(x, db, basis_t) {
    set_provenance(x, list(db = db, basis_t = basis_t,
                           timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

set_provenance <- function(x, prov) {
    if (is.null(prov)) return(x)
    attr(x, "patternq_provenance") <- prov
    x
}

# Put known columns first, in the given order, keeping any others after.
order_columns <- function(df, first) {
    first <- intersect(first, names(df))
    keep_provenance(df[, c(first, setdiff(names(df), first)), drop = FALSE], df)
}

#' Carry provenance over to a derived result
#'
#' Provenance is a plain attribute on results. Base R operations such as
#' \code{merge()} and \code{[} subsetting drop attributes, so provenance does
#' not follow a derived table automatically. patternq's own joins
#' (\code{add_*_context}) keep it; use this after your own merges or subsets.
#'
#' @param x Derived result
#' @param from The patternq result it came from
#' @return \code{x} with the provenance of \code{from}
#' @examples
#' \dontrun{
#' ctx <- keep_provenance(merge(measurements_df, other_df, by = "subject_id"), measurements_df)
#' }
#' @export
keep_provenance <- function(x, from) {
    set_provenance(x, provenance(from))
}

# Collapse cardinality-many values (list column) to "; "-joined strings.
join_many <- function(x, sep = "; ") {
    if (!is.list(x)) return(x)
    vapply(x, function(v) if (length(v) == 0 || all(is.na(v))) NA_character_
                          else paste(unlist(v), collapse = sep), character(1))
}

is_eid_map <- function(x) is_map(x) && identical(names(x), ":db/id")

# Pulled refs without a nested pattern come back as {":db/id": n}. Resolve
# those that are enums (have a :db/ident) to {":db/ident": ...} with one
# extra query, so enum values read as names whichever pull pattern was used.
resolve_enum_refs <- function(rows, db) {
    eids <- c()
    collect <- function(x) {
        if (is_eid_map(x)) { eids <<- c(eids, x[[":db/id"]]); return(invisible()) }
        if (is.list(x)) for (e in x) collect(e)
    }
    collect(rows)
    eids <- unique(eids)
    if (!length(eids)) return(rows)
    res <- raw_query(dq(find = c("?e", "?ident"), `in` = list(c("?e", "...")), args = list(I(eids)),
                        where = list(c("?e", ":db/ident", "?ident"))), db = db)
    if (!length(res$query_result)) return(rows)
    idents <- stats::setNames(vapply(res$query_result, function(r) as.character(r[[2]]), character(1)),
                              vapply(res$query_result, function(r) format(r[[1]], scientific = FALSE), character(1)))
    replace <- function(x) {
        if (is_eid_map(x)) {
            id <- idents[format(x[[":db/id"]], scientific = FALSE)]
            return(if (is.na(id)) x else list(":db/ident" = unname(id)))
        }
        if (is.list(x)) {
            nms <- names(x)
            x <- lapply(x, replace)
            names(x) <- nms
        }
        x
    }
    replace(rows)
}
