# API back-pressure contract (unify-central dev-docs/API_BACKPRESSURE.md),
# against mocked responses.

problem <- "urn:pattern-data-commons:problem:"

ok_resp <- function(headers = list(RateLimit = '"api";r=99;t=0')) {
    httr2::response(200, headers = c(list(`Content-Type` = "application/json"), headers),
                    body = charToRaw('{"query_result": [[1]], "basis_t": 7}'))
}

throttled <- function(kind = "rate-limited", status = 429, retry.after = NULL, retryable = TRUE) {
    headers <- list(`Content-Type` = "application/problem+json")
    if (!is.null(retry.after)) headers$`Retry-After` <- as.character(retry.after)
    body <- jsonlite::toJSON(list(type = paste0(problem, kind), title = "throttled", status = status,
                                  detail = paste(kind, "detail"), retryable = retryable),
                             auto_unbox = TRUE)
    httr2::response(status, headers = headers, body = charToRaw(body))
}

# Serve responses in order, recording the paths asked for and the sleeps.
local_stub <- function(responses, env = parent.frame()) {
    rec <- new.env()
    rec$paths <- character()
    rec$sleeps <- numeric()
    i <- 0
    httr2::local_mocked_responses(function(req) {
        rec$paths <- c(rec$paths, req$url)
        i <<- min(i + 1, length(responses))
        responses[[i]]
    }, env = env)
    local_mocked_bindings(pq_sleep = function(s) rec$sleeps <- c(rec$sleeps, s), .env = env)
    old.policy <- retry_policy()
    old.server <- set_query_server("http://stub.invalid")
    old.token <- set_token("test-token")
    withr::defer({
        do.call(set_retry_policy, old.policy)
        set_query_server(old.server)
        set_token(old.token)
        bp.env$next_allowed <- 0
        bp.env$max_timeout_ms <- default_max_timeout_ms
    }, envir = env)
    rec
}

q <- dq(find = "?x", where = list(c("?x", ":db/ident")))

test_that("429 and 503 are retried after Retry-After, with a message", {
    rec <- local_stub(list(throttled(retry.after = 3),
                           throttled("overloaded", status = 503, retry.after = 1), ok_resp()))
    msgs <- capture_messages(res <- raw_query(q, db = "db1", cache = FALSE))
    expect_match(msgs[1], "rate-limited")
    expect_match(msgs[2], "overloaded")
    expect_equal(res$basis_t, 7)
    expect_equal(rec$sleeps, c(3, 1))
    expect_length(rec$paths, 3)
})

test_that("without Retry-After the backoff is jittered and capped", {
    plain429 <- httr2::response(429, headers = list(`Content-Type` = "text/plain"), body = charToRaw("slow down"))
    rec <- local_stub(list(plain429, plain429, plain429, ok_resp()))
    set_retry_policy(max.backoff = 3)
    suppressMessages(raw_query(q, db = "db1", cache = FALSE))
    expect_length(rec$sleeps, 3)
    expect_true(all(rec$sleeps >= 0 & rec$sleeps <= pmin(3, 2^(1:3))))
})

test_that("gives up after max.retries with an error naming the problem", {
    rec <- local_stub(list(throttled(retry.after = 1)))
    set_retry_policy(max.retries = 2)
    err <- expect_error(suppressMessages(raw_query(q, db = "db1", cache = FALSE)),
                        class = "patternq_throttled")
    expect_equal(err$problem_type, paste0(problem, "rate-limited"))
    expect_match(conditionMessage(err), "rate-limited detail")
    expect_match(conditionMessage(err), "2 retries")
    expect_length(rec$paths, 3)
})

test_that("non-retryable problems, 502, 504 and 400 are not retried", {
    bad <- list(throttled(retryable = FALSE),
                httr2::response(502, body = charToRaw("bad gateway")),
                httr2::response(504, body = charToRaw("timeout")),
                httr2::response(400, headers = list(`Content-Type` = "application/json"),
                                body = charToRaw('{"error": "bad query"}')))
    for (r in bad) {
        rec <- local_stub(list(r, ok_resp()))
        expect_error(raw_query(q, db = "db1", cache = FALSE), "failed \\(HTTP")
        expect_length(rec$paths, 1)
        expect_length(rec$sleeps, 0)
    }
})

test_that("a query timeout says to narrow or page the query", {
    local_stub(list(httr2::response(400, headers = list(`Content-Type` = "application/json"),
                                    body = charToRaw(paste0('{"error": "Query canceled: timeout elapsed", ',
                                                            '"timeout": true, "type": "', problem, 'query-timeout"}')))))
    expect_error(raw_query(q, db = "db1", cache = FALSE), "Narrow the query or page it")
})

test_that("a query-too-broad refusal is surfaced and not retried", {
    rec <- local_stub(list(httr2::response(400, headers = list(`Content-Type` = "application/json"),
                                           body = charToRaw(paste0('{"error": "Query too broad: clause [?e ?a ?v] binds neither attribute nor entity; bind the attribute", "type": "', problem,
                                                                   'query-too-broad", "reason": "full-scan", "clause": "[?e ?a ?v]"}'))),
                           ok_resp()))
    expect_error(raw_query(q, db = "db1", cache = FALSE), "bind the attribute")
    expect_length(rec$paths, 1)
})

test_that("RateLimit r=0 makes the next call wait", {
    rec <- local_stub(list(ok_resp(list(RateLimit = '"api";r=0;t=1')), ok_resp()))
    raw_query(q, db = "db1", cache = FALSE)
    expect_length(rec$sleeps, 0)
    expect_message(raw_query(q, db = "db1", cache = FALSE), "rate limit reached")
    expect_length(rec$sleeps, 1)
    expect_true(rec$sleeps[1] > 0 && rec$sleeps[1] <= 1)
})

test_that("listing datasets and datoms go through back pressure", {
    rec <- local_stub(list(throttled(retry.after = 1),
                           httr2::response(200, headers = list(`Content-Type` = "application/json"),
                                           body = charToRaw('{"datasets": []}')),
                           throttled(retry.after = 1),
                           httr2::response(200, headers = list(`Content-Type` = "application/json"),
                                           body = charToRaw('{"datoms_chunk": [], "basis_t": 1}'))))
    suppressMessages({
        expect_equal(nrow(list_datasets()), 0)
        expect_equal(nrow(datoms("eavt", db = "db1")), 0)
    })
    expect_equal(rec$sleeps, c(1, 1))
})

test_that("timeouts above the advertised cap warn", {
    local_stub(list(ok_resp(list(`PDC-Query-Max-Timeout-Ms` = "60000"))))
    expect_warning(raw_query(q, db = "db1", cache = FALSE, timeout = 200), "cap of 120 s")
    expect_warning(raw_query(q, db = "db1", cache = FALSE, timeout = 90), "cap of 60 s")
})

test_that("policy defaults", {
    expect_equal(retry_policy(), list(max.retries = 5, max.backoff = 60, max.concurrency = 4))
})
