# Client side of the commons API back-pressure contract (unify-central
# dev-docs/API_BACKPRESSURE.md), the same in the R, Python, Clojure and Julia
# libraries:
#
# - 429 and 503 are retried, up to max.retries times, when the body is a
#   retryable application/problem+json, or is not problem+json at all (e.g.
#   from a load balancer). Every other status, 502 and 504 included, is
#   returned as is.
# - The wait is Retry-After when the server sends it, otherwise exponential
#   backoff with full jitter: runif(0, min(max.backoff, 2^attempt)) seconds.
# - Self-throttling: a RateLimit header with r=0 makes the next call wait t
#   seconds. An R session sends one request at a time, so it never exceeds
#   the server's per-key cap of 4 concurrent calls on its own.
# - Retries are reported with message(), with the problem type and the wait.
#
# Presigned S3 downloads are not API calls and don't go through here.

problem_prefix <- "urn:pattern-data-commons:problem:"
default_max_timeout_ms <- 120000

bp.env <- new.env(parent = emptyenv())
bp.env$retry_policy <- list(max.retries = 5, max.backoff = 60, max.concurrency = 4)
bp.env$next_allowed <- 0
bp.env$max_timeout_ms <- default_max_timeout_ms

# indirection so tests don't really sleep
pq_sleep <- function(seconds) Sys.sleep(seconds)

#' Set how this session handles API back pressure
#'
#' The commons API throttles heavy use: 429 or 503 responses with a
#' \code{Retry-After} header. patternq retries those calls, waiting
#' \code{Retry-After} (or a jittered exponential backoff when the server sends
#' none), and waits out the rate limit when the \code{RateLimit} header says
#' the key has no requests left. Retries are reported with \code{message()}.
#'
#' @param max.retries Retries per call before giving up with an error naming
#'   the problem type (default 5)
#' @param max.backoff Cap in seconds on the jittered backoff used when the
#'   server sends no \code{Retry-After} (default 60)
#' @param max.concurrency Concurrent \code{/query} and \code{/datoms} calls
#'   (default 4, the server's per-key cap). An R session sends one request at
#'   a time, so this only applies to parallel helpers.
#' @return The previous policy, invisibly
#' @export
set_retry_policy <- function(max.retries = NULL, max.backoff = NULL, max.concurrency = NULL) {
    old <- bp.env$retry_policy
    new <- old
    if (!is.null(max.retries)) new$max.retries <- as.integer(max.retries)
    if (!is.null(max.backoff)) new$max.backoff <- as.numeric(max.backoff)
    if (!is.null(max.concurrency)) {
        if (max.concurrency < 1) stop("max.concurrency must be at least 1", call. = FALSE)
        new$max.concurrency <- as.integer(max.concurrency)
    }
    bp.env$retry_policy <- new
    invisible(old)
}

#' The current API back-pressure policy
#'
#' @return A list with \code{max.retries}, \code{max.backoff} and
#'   \code{max.concurrency} (see \code{\link{set_retry_policy}})
#' @export
retry_policy <- function() bp.env$retry_policy

# The server's cap on requested query timeouts, from the last
# PDC-Query-Max-Timeout-Ms header seen (120000 until one is).
max_timeout_ms <- function() bp.env$max_timeout_ms

resp_problem <- function(resp) {
    if (!isTRUE(startsWith(httr2::resp_header(resp, "Content-Type", ""), "application/problem+json")))
        return(NULL)
    body <- tryCatch(jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE),
                     error = function(e) NULL)
    if (is.list(body)) body else NULL
}

resp_retry_after <- function(resp, problem) {
    h <- suppressWarnings(as.numeric(httr2::resp_header(resp, "Retry-After", NA)))
    if (!is.na(h)) return(max(0, h))
    if (is.numeric(problem$retry_after)) return(max(0, problem$retry_after))
    NULL
}

note_rate_headers <- function(resp) {
    rl <- httr2::resp_header(resp, "RateLimit")
    if (!is.null(rl)) {
        r <- regmatches(rl, regexec("\\br=([0-9]+)", rl))[[1]]
        t <- regmatches(rl, regexec("\\bt=([0-9]+)", rl))[[1]]
        if (length(r) == 2 && length(t) == 2 && as.numeric(r[2]) == 0)
            bp.env$next_allowed <- max(bp.env$next_allowed, now_seconds() + as.numeric(t[2]))
    }
    mt <- suppressWarnings(as.numeric(httr2::resp_header(resp, "PDC-Query-Max-Timeout-Ms", NA)))
    if (!is.na(mt)) bp.env$max_timeout_ms <- mt
}

now_seconds <- function() as.numeric(Sys.time())

wait_for_rate_limit <- function() {
    wait <- bp.env$next_allowed - now_seconds()
    if (wait > 0) {
        message(sprintf("patternq: rate limit reached (RateLimit r=0); waiting %.1f s before the next call", wait))
        pq_sleep(wait)
    }
}

# Perform an API request under the back-pressure contract. Returns the final
# response, which may still be an error for stop_for_response; stops with a
# patternq_throttled condition when throttling outlasts max.retries.
pq_perform <- function(req, what) {
    policy <- bp.env$retry_policy
    attempt <- 0
    repeat {
        wait_for_rate_limit()
        resp <- httr2::req_perform(req)
        note_rate_headers(resp)
        status <- httr2::resp_status(resp)
        if (!status %in% c(429, 503)) return(resp)
        problem <- resp_problem(resp)
        if (!is.null(problem) && !isTRUE(problem$retryable)) return(resp)
        ptype <- if (is.character(problem$type)) problem$type else NULL
        label <- if (is.null(ptype)) paste("HTTP", status) else ptype
        if (attempt >= policy$max.retries) {
            detail <- if (is.character(problem$detail)) problem$detail else
                substr(trimws(tryCatch(httr2::resp_body_string(resp), error = function(e) "")), 1, 500)
            rlang::abort(sprintf("%s throttled by the commons API (%s) and still throttled after %d retries: %s",
                                 what, label, attempt, detail),
                         class = "patternq_throttled", problem_type = ptype, status = status,
                         call = NULL)
        }
        attempt <- attempt + 1
        wait <- resp_retry_after(resp, problem)
        if (is.null(wait)) wait <- stats::runif(1, 0, min(policy$max.backoff, 2^attempt))
        message(sprintf("patternq: %s throttled (%s, HTTP %s); retry %d of %d in %.1f s",
                        what, label, status, attempt, policy$max.retries, wait))
        pq_sleep(wait)
    }
}
