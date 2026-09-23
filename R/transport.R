# HTTP transport to the Pattern Data Commons query service.
#
# Endpoints (all authenticated with a bearer API token):
#   POST /query/<db>           -> Accept text/plain: presigned URL of the gzipped, S3-cached
#                                 result {"query_result", "basis_t"}; Accept application/json:
#                                 the same JSON inline, skipping the cache
#   POST /datoms/<db>          -> JSON {"datoms_chunk": [...], "basis_t": ...}
#   POST /matrix/<db>/<key>    -> presigned URL of a gzipped TSV measurement matrix
#   GET  /api-v1/list/datasets -> JSON {"datasets": [...]}

pq_request <- function(path, accept = "text/plain") {
    httr2::request(query_server()) |>
        httr2::req_url_path_append(path) |>
        httr2::req_headers(Authorization = paste("Bearer", api_token()),
                           Accept = accept) |>
        httr2::req_user_agent("patternq-R") |>
        httr2::req_error(is_error = function(resp) FALSE)
}

stop_for_response <- function(resp, what) {
    status <- httr2::resp_status(resp)
    if (status == 200) return(invisible(NULL))
    body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    msg <- body
    parsed <- tryCatch(jsonlite::fromJSON(body), error = function(e) NULL)
    if (is.list(parsed) && !is.null(parsed$error)) msg <- parsed$error
    if (status == 401)
        msg <- paste("not authorized; check PATTERNQ_API_KEY / set_token().", msg)
    if (status == 403)
        msg <- paste("forbidden; your API key may not have access to this dataset.", msg)
    stop(sprintf("%s failed (HTTP %s): %s", what, status, msg), call. = FALSE)
}

# Download from a presigned S3 URL, gunzipping if needed.
fetch_presigned <- function(url) {
    resp <- httr2::request(trimws(url)) |>
        httr2::req_retry(max_tries = 5) |>
        httr2::req_error(is_error = function(resp) FALSE) |>
        httr2::req_perform()
    stop_for_response(resp, "Result download")
    raw <- httr2::resp_body_raw(resp)
    if (length(raw) >= 2 && raw[1] == as.raw(0x1f) && raw[2] == as.raw(0x8b))
        raw <- memDecompress(raw, type = "gzip")
    raw
}

query_body <- function(q, timeout, refresh.cache = FALSE) {
    qq <- q$query
    wire <- list()
    if (!is.null(qq$find)) wire[[":find"]] <- lapply(qq$find, wire_elem)
    ins <- c(list("$"), lapply(qq[["in"]], wire_elem))
    wire[[":in"]] <- ins
    for (sec in setdiff(names(qq), c("find", "in"))) {
        wire[[paste0(":", sec)]] <- lapply(qq[[sec]], wire_elem)
    }
    body <- list(query = wire, timeout = as.integer(timeout * 1000))
    if (refresh.cache) body[["refresh-cache"]] <- TRUE
    if (length(q$args)) body$args <- lapply(q$args, wire_arg)
    body
}

# Sequences (clauses, pull patterns, bindings) always serialize as JSON arrays;
# named lists serialize as JSON objects (pull pattern maps).
wire_seq <- function(x) {
    if (is.list(x) && !is.null(names(x)) && all(nzchar(names(x))))
        return(lapply(x, wire_seq))
    if (is.list(x))
        return(lapply(unname(x), wire_elem))
    as.list(unname(x))
}

wire_elem <- function(x) {
    if (is.atomic(x) && length(x) == 1 && !inherits(x, "AsIs")) return(unname(x))
    wire_seq(x)
}

wire_arg <- function(x) {
    if (is.data.frame(x))
        return(unname(lapply(seq_len(nrow(x)), function(i) unname(as.list(x[i, , drop = FALSE])))))
    if (is.atomic(x) && length(x) == 1 && !inherits(x, "AsIs")) return(unname(x))
    if (is.list(x)) return(lapply(unname(x), wire_arg))
    as.list(unname(unclass(x)))
}

#' Serialize a query to the JSON sent to the query service
#'
#' @param q A \code{patternq_query}
#' @param pretty Pretty-print
#' @param timeout Query timeout in seconds
#' @param refresh.cache Ask the service to recompute and re-cache the result
#' @return JSON string
#' @export
query_json <- function(q, pretty = FALSE, timeout = 30, refresh.cache = FALSE) {
    jsonlite::toJSON(query_body(q, timeout, refresh.cache), auto_unbox = TRUE, null = "null",
                     digits = NA, pretty = pretty)
}

#' Run a raw query and return the parsed response
#'
#' Low-level: returns the parsed JSON response (\code{query_result} and
#' \code{basis_t}) without converting to a data frame. Most users want
#' \code{\link{do_query}}.
#'
#' @param q A \code{patternq_query} (see \code{\link{dq}} and \code{\link{query}})
#' @param db Database name; defaults to \code{\link{current_db}}
#' @param timeout Query timeout in seconds
#' @param print.json Print the JSON request body
#' @param cache Use the service's S3 result cache (default). The service
#'   returns a presigned URL to a gzipped cached result, computing and caching
#'   it on a miss. With \code{FALSE} the result is returned inline as JSON and
#'   the cache is skipped entirely.
#' @param refresh.cache Recompute the result and overwrite the cached copy
#' @return A list with \code{query_result}, \code{basis_t} and \code{db}
#' @export
raw_query <- function(q, db = NULL, timeout = 30, print.json = FALSE,
                      cache = TRUE, refresh.cache = FALSE) {
    db <- ensure_db(db)
    q <- as_query(q)
    body <- query_json(q, timeout = timeout, refresh.cache = refresh.cache)
    if (print.json) message(body)
    accept <- if (cache) "text/plain" else "application/json"
    resp <- pq_request(c("query", db), accept = accept) |>
        httr2::req_body_raw(body, type = "application/json") |>
        httr2::req_timeout(timeout + 30) |>
        httr2::req_perform()
    stop_for_response(resp, "Query")
    payload <- httr2::resp_body_string(resp)
    if (startsWith(trimws(payload), "{")) {
        # errors (and some small results) come back inline
        res <- jsonlite::fromJSON(payload, simplifyVector = FALSE)
        if (!is.null(res$error)) stop("Query error: ", res$error, call. = FALSE)
    } else {
        raw <- fetch_presigned(payload)
        res <- jsonlite::fromJSON(rawToChar(raw), simplifyVector = FALSE)
    }
    if (!is.null(res$error)) stop("Query error: ", res$error, call. = FALSE)
    res$db <- db
    res
}

#' Run a query and return a data frame
#'
#' Relation results become a \code{data.frame} whose column names come from
#' the \code{:find} variables (\code{?sample-id} becomes \code{sample_id}).
#' Queries whose find spec is a single pull expression are flattened into
#' one column per attribute (see \code{\link{flatten_pull}}). The result
#' carries provenance (\code{\link{provenance}}): database name, basis t
#' and time of the query.
#'
#' @inheritParams raw_query
#' @param simplify If \code{TRUE}, a single-column relation result is returned as a vector
#' @param exclude.ids For pull results, drop \code{db/id} and \code{*/uid} columns
#' @return A \code{data.frame} (or vector, see \code{simplify})
#' @export
do_query <- function(q, db = NULL, timeout = 30, simplify = FALSE,
                     exclude.ids = TRUE, print.json = FALSE,
                     cache = TRUE, refresh.cache = FALSE) {
    q <- as_query(q)
    res <- raw_query(q, db = db, timeout = timeout, print.json = print.json,
                     cache = cache, refresh.cache = refresh.cache)
    rows <- res$query_result
    if (is_pull_query(q)) {
        rows <- resolve_enum_refs(rows, res$db)
        df <- pull_rows_to_df(rows, q, exclude.ids = exclude.ids)
    } else {
        df <- relation_to_df(rows, find_column_names(q))
    }
    df <- with_provenance(df, res$db, res$basis_t)
    if (simplify && is.data.frame(df) && ncol(df) == 1)
        return(with_provenance(df[[1]], res$db, res$basis_t))
    df
}

#' Run the same query against several databases and row-bind the results
#'
#' Each dataset lives in its own database, so cohort-level comparisons run the
#' same query per database. A \code{db} column is added to identify the source.
#'
#' @param q A query
#' @param dbs Character vector of database names
#' @param ... Passed to \code{\link{do_query}}
#' @return A \code{data.frame}
#' @export
across_dbs <- function(q, dbs, ...) {
    parts <- lapply(dbs, function(db) {
        df <- do_query(q, db = db, ...)
        if (nrow(df)) df$db <- db
        df
    })
    parts <- parts[vapply(parts, nrow, integer(1)) > 0]
    if (!length(parts)) return(data.frame())
    bind_rows_fill(parts)
}

#' Read raw datoms from an index
#'
#' @param index One of \code{"eavt"}, \code{"aevt"}, \code{"avet"}, \code{"vaet"}
#' @param components Index components, e.g. \code{list(":sample/id")}
#' @param db Database name
#' @param offset,limit Paging
#' @return A \code{data.frame} with columns \code{e}, \code{a}, \code{v}, \code{tx}
#' @export
datoms <- function(index, components = list(), db = NULL, offset = 0, limit = 1000) {
    db <- ensure_db(db)
    body <- list(index = paste0(":", sub("^:", "", index)),
                 components = as.list(components),
                 offset = offset, limit = limit)
    resp <- pq_request(c("datoms", db), accept = "application/json") |>
        httr2::req_body_json(body, auto_unbox = TRUE) |>
        httr2::req_perform()
    stop_for_response(resp, "Datoms request")
    res <- jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
    rows <- res$datoms_chunk
    df <- relation_to_df(lapply(rows, function(r) list(r[[":e"]], r[[":a"]], r[[":v"]], r[[":tx"]])),
                         c("e", "a", "v", "tx"))
    with_provenance(df, db, res$basis_t)
}

#' List datasets available to your API key
#'
#' @return A \code{data.frame} with one row per dataset: \code{dataset}
#'   (dataset name), \code{db} (current database name), \code{patient_count},
#'   \code{sample_count}, \code{assays} and \code{tags} (list columns)
#' @export
list_datasets <- function() {
    resp <- pq_request(c("api-v1", "list", "datasets"), accept = "application/json") |>
        httr2::req_perform()
    stop_for_response(resp, "Listing datasets")
    res <- jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
    ds <- res$datasets
    getv <- function(x, k, default = NA) if (is.null(x[[k]])) default else x[[k]]
    data.frame(
        dataset = vapply(ds, function(x) getv(x, "dataset/name", NA_character_), character(1)),
        db = vapply(ds, function(x) {
            d <- x[["dataset/database"]]
            if (is.null(d)) NA_character_ else getv(d, "database/name", NA_character_)
        }, character(1)),
        patient_count = vapply(ds, function(x) as.numeric(getv(x, "dataset/patient-count")), numeric(1)),
        sample_count = vapply(ds, function(x) as.numeric(getv(x, "dataset/sample-count")), numeric(1)),
        assays = I(lapply(ds, function(x) unlist(getv(x, "dataset/assays", list())))),
        tags = I(lapply(ds, function(x) unlist(getv(x, "dataset/tags", list())))),
        stringsAsFactors = FALSE
    )
}

#' Resolve a dataset name to its current database name
#'
#' Database names carry import dates and change when a dataset is re-imported;
#' dataset names are stable.
#'
#' @param dataset Dataset name, e.g. \code{"tcga-uvm"}
#' @return The database name
#' @export
resolve_db <- function(dataset) {
    ds <- list_datasets()
    hit <- ds$db[ds$dataset == dataset]
    if (!length(hit) || is.na(hit[1])) {
        # allow passing a database name straight through
        if (dataset %in% ds$db) return(dataset)
        stop(sprintf("Unknown dataset '%s'", dataset), call. = FALSE)
    }
    hit[1]
}

#' Download a measurement matrix
#'
#' Measurement matrices (e.g. single-cell counts) are stored as files rather
#' than as individual measurements. Find the keys with
#' \code{\link{measurement_matrices}}.
#'
#' @param matrix.key The matrix backing-file key
#' @param db Database name
#' @return A \code{data.frame} read from the (TSV) matrix file
#' @export
measurement_matrix <- function(matrix.key, db = NULL) {
    db <- ensure_db(db)
    resp <- pq_request(c("matrix", db, matrix.key)) |>
        httr2::req_body_json(stats::setNames(list(), character(0))) |>
        httr2::req_perform()
    stop_for_response(resp, "Matrix request")
    raw <- fetch_presigned(httr2::resp_body_string(resp))
    utils::read.delim(text = rawToChar(raw), check.names = FALSE, stringsAsFactors = FALSE)
}
