# The transit formats are optional for users but always tested: a missing
# transit package fails here instead of skipping.
stopifnot("the transit tests need transit-r: remotes::install_github(\"vendekagon-labs/transit-r\")" =
              requireNamespace("transit", quietly = TRUE))

test_that("transit responses decode to the JSON response's shape", {
    json <- paste0('{"query_result":[[{":sample/id":"S1",":sample/type":{":db/ident":":sample.type/tumor"},',
                   '":sample/recurrence":false,":sample/uid":["A","B"]},true,17592186270456,1.5,null,',
                   '"2026-03-31T22:21:35.17Z","2026-03-31T22:25:15.903Z","1970-01-01T00:00:00Z"]],"basis_t":42}')
    tj <- paste0('["^ ","query_result",[[["^ ","~:sample/id","S1","~:sample/type",["^ ","~:db/ident","~:sample.type/tumor"],',
                 '"~:sample/recurrence",false,"~:sample/uid",["A","B"]],true,17592186270456,1.5,null,',
                 '"~m1774995695170","~m1774995915903","~m0"]],"basis_t",42]')
    expected <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    expect_identical(decode_transit(charToRaw(tj), "transit+json"), expected)
    value <- transit::from_transit(tj)
    mp <- transit::to_transit(value, "msgpack")
    expect_identical(decode_transit(mp, "transit+msgpack"), expected)
    expect_error(check_format("edn"), "format must be one of")
})

test_that("transit formats give the same results as JSON", {
    skip_if_no_token()
    db <- test_db("tcga-uvm")
    # pulled attributes may come in a different column order (map key order)
    strip <- function(df) {
        attr(df, "patternq_provenance") <- NULL
        df <- df[, sort(names(df))]
        df <- df[do.call(order, unname(as.list(df[vapply(df, is.atomic, NA)]))), ]
        rownames(df) <- NULL
        df
    }
    for (f in list(samples, variants)) {
        a <- strip(f(db, cache = FALSE))
        expect_identical(strip(f(db, format = "transit+json")), a)
        expect_equal(strip(f(db, format = "transit+msgpack")), a, tolerance = 1e-6)
    }
})
