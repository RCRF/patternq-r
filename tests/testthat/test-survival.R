test_that("logrank_test matches survival::survdiff", {
    skip_if_not_installed("survival")
    set.seed(1)
    t <- rexp(60, 0.1); e <- runif(60) > 0.3; g <- sample(c("a", "b", "c"), 60, TRUE)
    lr <- logrank_test(t, e, g)
    sd <- survival::survdiff(survival::Surv(t, e) ~ g)
    expect_equal(lr$chisq, sd$chisq, tolerance = 1e-10)
    expect_equal(lr$df, 2)
})

test_that("logrank_test fixed reference (shared with Python and Clojure tests)", {
    time <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    event <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
    group <- c("a", "b", "a", "a", "b", "a", "b", "b", "a", "b")
    expect_equal(round(logrank_test(time, event, group)$chisq, 6), round(0.2201257417, 6))
})

test_that("median_split and survival_status", {
    expect_equal(median_split(c(1, 2, 3, 4, NA)), c("low", "low", "high", "high", NA))
    oc <- data.frame(os = c(20, 5, 5, NA), os_event = c(FALSE, TRUE, FALSE, TRUE))
    expect_equal(survival_status(oc), c("alive at 1 year", "died within 1 year", NA, NA))
})

test_that("change_from_baseline", {
    tab <- data.frame(subject_id = c("p1", "p1", "p2", "p2", "p3"),
                      timepoint_id = c("C1D1", "C2D1", "C1D1", "C2D1", "C2D1"),
                      cell_population = "T", value = c(1, 4, 2, 2, 9))
    ch <- change_from_baseline(tab)
    expect_equal(nrow(ch), 4)  # p3 has no baseline
    expect_equal(ch$change[ch$subject_id == "p1" & ch$timepoint_id == "C2D1"], 2)
    expect_equal(with(change_from_baseline(tab, method = "difference"), change[subject_id == "p1" & timepoint_id == "C2D1"]), 3)
})

test_that("survival_by_median attaches a log-rank test", {
    oc <- data.frame(subject_id = paste0("s", 1:8), os = c(1, 2, 3, 4, 10, 11, 12, 13),
                     os_event = TRUE)
    vals <- setNames(c(1, 2, 3, 4, 5, 6, 7, 8), paste0("s", 1:8))
    km <- survival_by_median(vals, oc)
    expect_equal(sort(unique(km$group)), c("high", "low"))
    expect_lt(attr(km, "logrank")$p, 0.05)
})
