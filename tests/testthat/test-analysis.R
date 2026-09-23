test_that("percentile_rank uses mid-rank ties", {
    expect_equal(percentile_rank(1:4, 2.5), 50)
    expect_equal(percentile_rank(c(1, 2, 2, 3), 2), 50)
    expect_equal(percentile_rank(1:4, 10), 100)
    expect_true(is.na(percentile_rank(numeric(0), 1)))
})

test_that("log_fold_change matches the Clojure definition", {
    a <- c(G1 = 0, G2 = 3, G3 = 7)
    b <- c(G1 = 1, G2 = 3, G4 = 15)
    lfc <- log_fold_change(a, b)
    expect_equal(lfc$lfc[lfc$hgnc_symbol == "G1"], 1)
    expect_equal(lfc$lfc[lfc$hgnc_symbol == "G2"], 0)
    expect_equal(lfc$lfc[lfc$hgnc_symbol == "G3"], -3)
    expect_equal(lfc$lfc[lfc$hgnc_symbol == "G4"], 4)
    expect_equal(lfc$avg_log10[lfc$hgnc_symbol == "G2"], log10(4))
})

test_that("ssgsea_score: higher when the set is highly expressed", {
    x <- stats::setNames(seq(100, 1), paste0("G", 1:100))
    top <- ssgsea_score(x, paste0("G", 1:10))
    bottom <- ssgsea_score(x, paste0("G", 91:100))
    expect_gt(top, 0)
    expect_lt(bottom, 0)
    expect_true(is.na(ssgsea_score(x, "nope")))
})

test_that("expression distances", {
    a <- c(A = 1, B = 0); b <- c(A = 0, B = 1)
    expect_equal(expression_distance(a, a), 0)
    expect_equal(expression_distance(a, b), 1)
    expect_equal(expression_distance(a, b, "euclidean"), sqrt(2))
    m <- rbind(s1 = c(A = 1, B = 0), s2 = c(A = 0, B = 1))
    expect_equal(nearest_samples(a, m, n = 1)$sample_id, "s1")
})

test_that("top_varying_genes and kaplan_meier", {
    m <- cbind(flat = rep(5, 4), varied = c(0, 10, 100, 1000))
    expect_equal(top_varying_genes(m, 1)$hgnc_symbol, "varied")
    km <- kaplan_meier(c(1, 2, 3, 4), c(TRUE, FALSE, TRUE, TRUE))
    expect_equal(km$surv, c(0.75, 0.75, 0.375, 0))
})

test_that("genesets ship with the package", {
    gs <- genesets()
    expect_true(all(c("hallmark_hypoxia", "housekeeping", "multi_cancer") %in% names(gs)))
    expect_true("GAPDH" %in% geneset("housekeeping"))
})

test_that("live: sample vs cohort comparison and two-sample change", {
    skip_if_no_token()
    db <- test_db("H37001"); uvm <- test_db("tcga-uvm")
    genes <- c("MLANA", "PMEL", "TYR", "GAPDH", "BAP1", "PRAME")
    cmp <- compare_to_cohort("H37001-003", db = db, cohort.db = uvm, genes = genes)
    expect_setequal(cmp$hgnc_symbol, genes)
    expect_true(all(cmp$cohort_n == 80))
    expect_true(all(cmp$percentile >= 0 & cmp$percentile <= 100))
    expect_s3_class(plot_zscores(cmp), "plotly")
    ch <- compare_samples("H37001-003", "H37001-001", db = db)
    expect_true(all(ch$avg_log10 >= 0.5))
    expect_s3_class(plot_ma(ch), "plotly")
    expect_s3_class(examine_geneset(c("MLANA", "PMEL"), "H37001-003", db = db, cohort.dbs = uvm), "plotly")
})
