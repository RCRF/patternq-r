# Live tests query the dev query service; they run when PATTERNQ_API_KEY is set.
skip_if_no_token <- function() {
    if (!nzchar(Sys.getenv("PATTERNQ_API_KEY"))) skip("PATTERNQ_API_KEY not set")
}

test_db <- local({
    cache <- list()
    function(dataset) {
        if (is.null(cache[[dataset]])) cache[[dataset]] <<- resolve_db(dataset)
        cache[[dataset]]
    }
})
