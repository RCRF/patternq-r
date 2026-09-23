# Query objects.
#
# A patternq query is list(query = list(find=, in=, where=, ...), args = list())
# holding the same structure as the JSON sent to the query service:
# variables are "?x" strings, attributes/keywords are ":ns/name" strings,
# clauses are vectors/lists, pull patterns use named lists for maps.
# The implicit database "$" is always prepended to :in when the query runs.

#' Build a query from plain R data
#'
#' The query mirrors the JSON wire format one-to-one, so any Datomic query can
#' be written directly: variables are \code{"?x"} strings, attributes are
#' \code{":ns/name"} strings, each clause is a character vector (or a list
#' when it nests), and pull pattern maps are named lists.
#'
#' Arguments are bound in order to the \code{in} inputs. A length-one value
#' binds a scalar; a longer vector (or one wrapped in \code{I()}) binds a
#' collection input such as \code{list("?ids", "...")}.
#'
#' Service limitation: string literals starting with \code{:} or \code{?}
#' cannot be expressed in a query; pass them as arguments.
#'
#' @param find The find spec: a character vector of variables, or a list whose
#'   elements are variables, aggregates (\code{c("count", "?x")}) or pull
#'   expressions (\code{list("pull", "?e", c("*"))})
#' @param where List of clauses
#' @param in Inputs other than \code{$}, e.g. \code{list("?sample-id")} or
#'   \code{list(c("?gene", "..."))}
#' @param args Query arguments, bound to \code{in} in order
#' @param with Optional \code{:with} variables
#' @return A \code{patternq_query}
#' @examples
#' q <- dq(find = c("?id", "?sex"),
#'         where = list(c("?s", ":subject/id", "?id"),
#'                      c("?s", ":subject/sex", "?sx"),
#'                      c("?sx", ":db/ident", "?sex")))
#' @export
dq <- function(find, where, `in` = NULL, args = list(), with = NULL) {
    q <- list(find = as.list(find), `in` = if (is.null(`in`)) NULL else as.list(`in`),
              with = if (is.null(with)) NULL else as.list(with),
              where = as.list(where))
    new_query(q, args = args)
}

new_query <- function(q, args = list()) {
    q <- q[!vapply(q, is.null, logical(1))]
    structure(list(query = q, args = as.list(args)), class = "patternq_query")
}

#' Coerce to a patternq query
#'
#' Accepts a \code{patternq_query}, a query built with the datalogr-style DSL
#' (\code{\link{query}}), or a list with \code{find}/\code{where} (optionally
#' \code{in}, \code{args}) elements.
#'
#' @param q Query-like object
#' @return A \code{patternq_query}
#' @export
as_query <- function(q) {
    if (inherits(q, "patternq_query")) return(q)
    if (is.list(q) && !is.null(q$query)) {
        qq <- q$query
        names(qq) <- sub("^:", "", names(qq))
        args <- q$args
        if (is.null(args)) args <- list()
        return(new_query(qq, args = args))
    }
    if (is.list(q) && !is.null(q$find)) {
        names(q) <- sub("^:", "", names(q))
        args <- q$args
        q$args <- NULL
        if (is.null(args)) args <- list()
        return(new_query(q, args = args))
    }
    stop("Not a query: expected dq(...), query(...) or list(find = ..., where = ...)", call. = FALSE)
}

#' Add arguments to a query
#'
#' @param q A query
#' @param ... Arguments appended in order
#' @return The query with arguments appended
#' @export
with_args <- function(q, ...) {
    q <- as_query(q)
    q$args <- c(q$args, list(...))
    q
}

#' @export
print.patternq_query <- function(x, ...) {
    cat("<patternq query>\n")
    cat(query_json(x, pretty = TRUE), "\n")
    invisible(x)
}

#' Does a query's find spec contain a pull expression?
#' @param q A query
#' @return Logical
#' @export
is_pull_query <- function(q) {
    q <- as_query(q)
    any(vapply(q$query$find, function(e) {
        !is.atomic(e) || length(e) > 1
    }, logical(1)) & vapply(q$query$find, function(e) {
        identical(as.character(unlist(e)[1]), "pull")
    }, logical(1)))
}

#' Convert a Datalog variable or keyword to a column name
#'
#' \code{"?sample-id"} becomes \code{"sample_id"}; \code{":subject/id"}
#' becomes \code{"subject_id"}.
#'
#' @param x Character vector
#' @return Character vector of syntactic, snake_case names
#' @export
clean_names <- function(x) {
    x <- sub("^[?:]", "", x)
    x <- gsub("[-/.]", "_", x)
    make.names(x)
}

find_column_names <- function(q) {
    nms <- vapply(q$query$find, function(e) {
        if (is.atomic(e) && length(e) == 1) return(as.character(e))
        parts <- as.character(unlist(e))
        paste(sub("^[?:]", "", parts), collapse = "_")
    }, character(1))
    make.unique(clean_names(nms), sep = "_")
}
