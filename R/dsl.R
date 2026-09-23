# Datalog DSL, ported from datalogr (Parker Institute for Cancer
# Immunotherapy; relicensed Apache-2.0).
#
# Writes Datomic queries with R syntax:
#
#   query(
#     find(?sample-id, ?vaf),
#     where(
#       d(?m, measurement/vaf, ?vaf),
#       d(?m, measurement/sample, ?s),
#       d(?s, sample/id, ?sample-id)
#     ),
#     args(?gene <- "TP53")
#   )
#
# Conversion rules:
#   - ?name symbols are Datalog variables (?gene-hgnc-symbol is fine)
#   - ns/name is an attribute keyword (:ns/name)
#   - . is _ (blank)
#   - infix < > become prefix expression clauses
#   - f(x) becomes an expression (f x)
#   - !!x substitutes the value of the R variable x
#
# The query object produced is a patternq_query (see dq()), so DSL queries,
# dq() queries and canned queries all run through do_query().

process_question_mark <- function(s) {
    ret <- gsub("`", "", s)
    if (grepl("?", ret, fixed = TRUE)) {
        regex <- "\\?\\(([^()]*)\\)"
        m <- regmatches(ret, gregexpr(regex, ret))[[1]]
        for (x in m) {
            inner <- sub(regex, "\\1", x)
            inner <- gsub(" - ", "-", inner, fixed = TRUE)
            ret <- sub(x, paste0("?", inner), ret, fixed = TRUE)
        }
    }
    trimws(ret)
}

try_numeric <- function(s) {
    if (!is.character(s)) return(s)
    out <- suppressWarnings(as.numeric(s))
    if (length(s) == 1 && !is.na(out)) out else trimws(s)
}

parse_expression <- function(s) {
    ret <- s
    if (!is.character(ret)) {
        if (is.numeric(ret) || is.logical(ret)) return(ret)
        ret <- paste(deparse(ret), collapse = "")
        ret <- process_question_mark(ret)

        # keyword
        if (grepl("[a-z]/[a-z._]", ret, ignore.case = TRUE) && !grepl("^\\?", ret) &&
            !grepl("\\(", ret)) {
            ret <- gsub(" ", "", ret)
            ret <- paste0(":", ret)
            return(ret)
        }

        # infix comparison
        if (grepl("<|>", ret)) {
            op <- regmatches(ret, gregexpr("<=|>=|<|>", ret))[[1]][1]
            v <- trimws(strsplit(ret, op, fixed = TRUE)[[1]])
            return(list(list(op, try_numeric(v[[1]]), try_numeric(v[[2]]))))
        }

        # vector c(?a, ?b) -> list(?a, ?b)
        if (startsWith(ret, "c(")) {
            m <- gsub(" ", "", ret)
            inner <- sub("^c\\((.*)\\)$", "\\1", m)
            return(lapply(as.list(strsplit(inner, ",")[[1]]), try_numeric))
        }

        # function call f(x, y) -> (f x y)
        if (grepl("^[^(]+\\(.*\\)$", ret)) {
            f <- sub("^([^(]+)\\(.*\\)$", "\\1", ret)
            a <- sub("^[^(]+\\((.*)\\)$", "\\1", ret)
            a <- trimws(strsplit(a, ",")[[1]])
            return(lapply(c(f, a), try_numeric))
        }

        if (ret == ".") ret <- "_"
        return(try_numeric(ret))
    }
    # character literals (e.g. "001") are kept as strings
    ret
}

#' Convert R expressions to Datalog terms
#'
#' Builds one Datalog clause from R expressions (see the package DSL notes in
#' \code{\link{query}}).
#'
#' @param ... Expressions forming one clause
#' @return A list representing the clause
#' @examples
#' d(?m, measurement/sample, ?s)
#' @export
d <- function(...) {
    ret <- rlang::exprs(...)
    ret <- lapply(ret, parse_expression)
    names(ret) <- NULL
    # expression clauses come back wrapped: [[op a b]]
    if (length(ret) == 1 && is.list(ret[[1]]) && length(ret[[1]]) == 1 && is.list(ret[[1]][[1]]))
        return(ret[[1]])
    ret
}

#' The where portion of a query
#' @param ... Clauses built with \code{\link{d}}, \code{\link{or}},
#'   \code{\link{not}}, \code{\link{not_join}}
#' @return Query section
#' @export
where <- function(...) {
    ret <- lapply(rlang::list2(...), identity)
    names(ret) <- NULL
    list(where = ret)
}

#' The find portion of a query
#' @param ... Variables, aggregates (e.g. \code{count(?s)}) or \code{\link{pull}}
#'   expressions
#' @return Query section
#' @export
find <- function(...) {
    exprs <- rlang::exprs(...)
    env <- parent.frame()
    ret <- lapply(exprs, function(e) {
        if (startsWith(paste(deparse(e), collapse = "")[[1]], "pull"))
            eval(e, env)
        else
            parse_expression(e)
    })
    names(ret) <- NULL
    list(find = ret)
}

#' Construct a Datalog query
#'
#' Combines \code{\link{find}}, \code{\link{where}} and optionally
#' \code{\link{args}} and \code{\link{with_vars}} sections into a
#' \code{patternq_query}.
#'
#' @param ... Query sections
#' @return A \code{patternq_query}
#' @export
query <- function(...) {
    parts <- list(...)
    q <- list()
    args <- list()
    for (p in parts) {
        for (k in names(p)) {
            if (k == "args") args <- c(args, p[[k]])
            else q[[k]] <- c(q[[k]], p[[k]])
        }
    }
    new_query(q, args = args)
}

#' Query inputs and their values
#'
#' Each argument has the form \code{?variable <- value}. The value is
#' evaluated in the calling environment. A value of length greater than one
#' (or wrapped in \code{I()}) is bound as a collection (\code{[?variable ...]}).
#'
#' @param ... Bindings
#' @return Query section
#' @examples
#' args(?sample-id <- "H37001-001", ?gene <- c("TP53", "KRAS"))
#' @export
args <- function(...) {
    exprs <- rlang::exprs(...)
    env <- parent.frame()
    ins <- list()
    vals <- list()
    for (e in exprs) {
        # `?` has the lowest precedence in R, so `?x <- v` parses as `?`(x <- v)
        if (is.call(e) && identical(e[[1]], as.name("?")) && length(e) == 2 &&
            is.call(e[[2]]) && identical(e[[2]][[1]], as.name("<-"))) {
            b <- e[[2]]
            var <- paste0("?", gsub(" - ", "-", paste(deparse(b[[2]]), collapse = ""), fixed = TRUE))
            val <- eval(b[[3]], env)
        } else if (is.call(e) && identical(e[[1]], as.name("<-"))) {
            var <- process_question_mark(paste(deparse(e[[2]]), collapse = ""))
            val <- eval(e[[3]], env)
        } else {
            stop("args() expects bindings of the form ?var <- value", call. = FALSE)
        }
        if (length(val) > 1 || inherits(val, "AsIs") || is.list(val)) {
            ins[[length(ins) + 1]] <- list(var, "...")
            vals[[length(vals) + 1]] <- I(unclass(val))
        } else {
            ins[[length(ins) + 1]] <- var
            vals[[length(vals) + 1]] <- val
        }
    }
    list(`in` = ins, args = vals)
}

#' or clause
#' @param ... Clauses built with \code{\link{d}}
#' @return Clause
#' @export
or <- function(...) {
    c(list("or"), unname(list(...)))
}

#' Generate an or clause over alternative values
#'
#' \code{generate_or(c("A", "B"), ?p, epitope/id)} is equivalent to
#' \code{or(d(?p, epitope/id, "A"), d(?p, epitope/id, "B"))}.
#'
#' @param alternatives Values
#' @param ... Partial clause
#' @return Clause
#' @export
generate_or <- function(alternatives, ...) {
    partial <- d(...)
    c(list("or"), lapply(alternatives, function(x) c(partial, list(x))))
}

#' not clause
#' @inheritParams or
#' @return Clause
#' @export
not <- function(...) {
    c(list("not"), unname(list(...)))
}

#' not-join clause
#' @param vars Variables to join on, built with \code{\link{d}}
#' @param ... Clauses
#' @return Clause
#' @export
not_join <- function(vars, ...) {
    c(list("not-join", vars), unname(list(...)))
}

#' with portion of a query
#'
#' Named \code{with_vars} (datalogr called it \code{with}) so that
#' \code{base::with} is not masked.
#'
#' @param ... Variables
#' @return Query section
#' @export
with_vars <- function(...) {
    list(with = d(...))
}

#' Combine queries
#'
#' Concatenates the sections and arguments of two queries.
#'
#' @param a,b Queries
#' @return A \code{patternq_query}
#' @export
c_query <- function(a, b) {
    a <- as_query(a)
    b <- as_query(b)
    for (f in union(names(a$query), names(b$query)))
        a$query[[f]] <- c(a$query[[f]], b$query[[f]])
    a$args <- c(a$args, b$args)
    a
}

#' Datomic pull expression
#'
#' Rules (as in datalogr): vectors are written \code{c()}, maps
#' \code{{key = val}}, the wildcard \code{*} is written \code{.}, reverse
#' references use \code{.} instead of \code{_} (\code{variant/.gene}).
#'
#' @param ... Exactly two arguments: the entity variable and the pattern
#' @return A find element
#' @examples
#' pull(?s, c(., {sample/subject = c(subject/id)}))
#' @export
pull <- function(...) {
    exprs <- rlang::exprs(...)
    if (length(exprs) != 2)
        stop("pull() takes exactly 2 arguments", call. = FALSE)
    v <- process_question_mark(paste(deparse(exprs[[1]]), collapse = ""))
    list("pull", v, parse_pull_pattern(exprs[[2]]))
}

parse_pull_pattern <- function(e) {
    if (is.call(e) && (identical(e[[1]], as.name("c")) || identical(e[[1]], as.name("{")))) {
        is_map <- identical(e[[1]], as.name("{"))
        items <- as.list(e)[-1]
        if (is_map) {
            out <- list()
            for (it in items) {
                # {a/b = c(...)} parses as `=`(a/b, c(...))
                key <- pull_attr(it[[2]])
                out[[key]] <- parse_pull_pattern(it[[3]])
            }
            return(out)
        }
        out <- lapply(items, function(it) {
            if (is.call(it) && identical(it[[1]], as.name("{"))) {
                return(parse_pull_pattern(it))
            }
            if (is.call(it) && identical(it[[1]], as.name("c")))
                return(parse_pull_pattern(it))
            pull_attr(it)
        })
        return(out)
    }
    list(pull_attr(e))
}

pull_attr <- function(e) {
    s <- gsub(" ", "", paste(deparse(e), collapse = ""))
    if (s == ".") return("*")
    if (s == "as") return(":as")
    if (grepl("/", s, fixed = TRUE)) {
        s <- sub("/\\.", "/_", s)
        return(paste0(":", s))
    }
    try_numeric(s)
}
