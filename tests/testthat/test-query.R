test_that("dq serializes to the service's JSON form", {
    q <- dq(find = c("?id"), where = list(c("?s", ":subject/id", "?id")))
    j <- jsonlite::fromJSON(query_json(q), simplifyVector = FALSE)
    expect_equal(j$query[[":find"]], list("?id"))
    expect_equal(j$query[[":in"]], list("$"))
    expect_equal(j$query[[":where"]], list(list("?s", ":subject/id", "?id")))
})

test_that("collection args bind as [?x ...]", {
    q <- dq(find = "?s", where = list(c("?s", ":sample/id", "?id")),
            `in` = list(c("?id", "...")), args = list(I("A")))
    j <- jsonlite::fromJSON(query_json(q), simplifyVector = FALSE)
    expect_equal(j$query[[":in"]], list("$", list("?id", "...")))
    expect_equal(j$args, list(list("A")))
})

test_that("DSL builds the same query as dq", {
    genes <- c("TP53", "KRAS")
    q1 <- query(find(?sample-id, ?vaf),
                where(d(?m, measurement/vaf, ?vaf),
                      d(?m, measurement/sample, ?s),
                      d(?s, sample/id, ?sample-id),
                      d(?vaf > 0.3)),
                args(?gene <- genes))
    j <- jsonlite::fromJSON(query_json(q1), simplifyVector = FALSE)
    expect_equal(j$query[[":find"]], list("?sample-id", "?vaf"))
    expect_equal(j$query[[":where"]][[4]], list(list(">", "?vaf", 0.3)))
    expect_equal(j$query[[":in"]], list("$", list("?gene", "...")))
    expect_equal(j$args, list(list("TP53", "KRAS")))
})

test_that("pull patterns", {
    p <- pull(?s, c(., {sample/subject = c(subject/id)}))
    expect_equal(p, list("pull", "?s", list("*", list(":sample/subject" = list(":subject/id")))))
})

test_that("flatten_pull names columns after attributes and resolves enums", {
    m <- list(":sample/id" = "S1", ":db/id" = 1,
              ":sample/subject" = list(":subject/id" = "P1"),
              ":sample/specimen" = list(":db/ident" = ":sample.specimen/ffpe"))
    f <- flatten_pull(m)
    expect_equal(f, list(sample_id = "S1", subject_id = "P1", sample_specimen = "ffpe"))
})

test_that("aggregate column names", {
    q <- dq(find = list(c("count", "?s")), where = list(c("?s", ":sample/id")))
    expect_equal(patternq:::find_column_names(q), "count_s")
})

test_that("keep_provenance carries provenance to derived tables", {
    df <- patternq:::with_provenance(data.frame(subject_id = c("a", "b"), v = 1:2), "db1", 42)
    expect_identical(class(df), "data.frame")
    m <- merge(df, data.frame(subject_id = "a", w = 3), by = "subject_id")
    expect_null(provenance(m))
    expect_equal(provenance(keep_provenance(m, df))$db, "db1")
})
