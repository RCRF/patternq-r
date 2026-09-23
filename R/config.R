#' patternq: query and analysis tools for the Pattern Data Commons
#'
#' Configure access with the \code{PATTERNQ_ENDPOINT} and
#' \code{PATTERNQ_API_KEY} environment variables, or in-session with
#' \code{\link{set_query_server}} and \code{\link{set_token}}. Each dataset
#' lives in its own database: select one with \code{\link{set_db}} or pass
#' \code{db} to any query function.
#'
#' @keywords internal
"_PACKAGE"

pkg.env <- new.env(parent = emptyenv())

default_endpoint <- "https://data-commons.rcrf-dev.org"

#' Query service endpoint in use
#'
#' Resolution order: value set with \code{\link{set_query_server}}, then the
#' \code{PATTERNQ_ENDPOINT} environment variable, then the default dev server.
#'
#' @return The endpoint URL as a string
#' @export
query_server <- function() {
    if (!is.null(pkg.env$endpoint))
        return(pkg.env$endpoint)
    env <- Sys.getenv("PATTERNQ_ENDPOINT", "")
    if (startsWith(env, "http"))
        return(sub("/+$", "", env))
    default_endpoint
}

#' Set the query service endpoint for this session
#'
#' @param url The endpoint URL, e.g. \code{"https://data-commons.rcrf-dev.org"}
#' @return The previous value, invisibly
#' @export
set_query_server <- function(url) {
    old <- pkg.env$endpoint
    pkg.env$endpoint <- if (is.null(url)) NULL else sub("/+$", "", url)
    invisible(old)
}

#' Set the API token for this session
#'
#' Overrides the \code{PATTERNQ_API_KEY} environment variable. Obtain a key
#' from the user settings page of the Pattern Data Commons dashboard.
#'
#' @param token The API token
#' @return The previous value, invisibly
#' @export
set_token <- function(token) {
    old <- pkg.env$token
    pkg.env$token <- token
    invisible(old)
}

api_token <- function() {
    token <- pkg.env$token
    if (is.null(token) || !nzchar(token))
        token <- Sys.getenv("PATTERNQ_API_KEY", "")
    if (!nzchar(token))
        stop("No API token: set PATTERNQ_API_KEY in the environment or call set_token()",
             call. = FALSE)
    token
}

#' Set the default database for this session
#'
#' Every dataset is its own database, so the database name selects the
#' dataset. Use \code{\link{resolve_db}} to go from a dataset name (e.g.
#' \code{"tcga-uvm"}) to its current database name.
#'
#' @param db Database name, e.g. \code{"tcga-gdc-uvm-3"}
#' @return The previous value, invisibly
#' @export
set_db <- function(db) {
    old <- pkg.env$db
    pkg.env$db <- db
    invisible(old)
}

#' The default database for this session, if set
#' @return The database name or \code{NULL}
#' @export
current_db <- function() {
    pkg.env$db
}

ensure_db <- function(db) {
    if (is.null(db))
        db <- pkg.env$db
    if (is.null(db))
        stop("No database given: pass db = ... or call set_db()", call. = FALSE)
    db
}
