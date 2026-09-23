# Survival helpers for outcome-association analyses (as in the PRINCE
# biomarker figures: OS stratified at the median of a baseline biomarker,
# log-rank p-values, landmark survival status).

#' Kaplan-Meier estimate
#' @param time Numeric times
#' @param event Logical, \code{TRUE} = event, \code{FALSE} = censored
#' @return \code{data.frame} of \code{time}, \code{n.risk}, \code{n.event},
#'   \code{n.censor}, \code{surv}
#' @export
kaplan_meier <- function(time, event) {
    event <- as.logical(event)
    o <- order(time)
    time <- time[o]; event <- event[o]
    ut <- unique(time)
    n.risk <- vapply(ut, function(t) sum(time >= t), numeric(1))
    n.event <- vapply(ut, function(t) sum(time == t & event), numeric(1))
    n.censor <- vapply(ut, function(t) sum(time == t & !event), numeric(1))
    surv <- cumprod(1 - n.event / n.risk)
    data.frame(time = ut, n.risk = n.risk, n.event = n.event, n.censor = n.censor, surv = surv)
}

#' Log-rank test
#'
#' Mantel-Haenszel log-rank test for a difference in survival between two or
#' more groups (hand-rolled; no survival package). Matches
#' \code{survival::survdiff}.
#'
#' @param time Numeric times
#' @param event Logical or 0/1 events
#' @param group Group labels
#' @return List: \code{chisq}, \code{df}, \code{p}, \code{observed},
#'   \code{expected} (per group), \code{n} (per group)
#' @export
logrank_test <- function(time, event, group) {
    keep <- !is.na(time) & !is.na(event) & !is.na(group)
    time <- time[keep]; event <- as.logical(event[keep]); group <- as.character(group[keep])
    lv <- sort(unique(group))
    k <- length(lv)
    if (k < 2) return(list(chisq = NA_real_, df = 0, p = NA_real_))
    times <- sort(unique(time[event]))
    O <- numeric(k); E <- numeric(k); V <- matrix(0, k, k)
    for (t in times) {
        at.risk <- vapply(lv, function(g) sum(time >= t & group == g), numeric(1))
        d.g <- vapply(lv, function(g) sum(time == t & event & group == g), numeric(1))
        n <- sum(at.risk); d <- sum(d.g)
        if (n < 1) next
        O <- O + d.g
        E <- E + d * at.risk / n
        if (n > 1) {
            f <- d * (n - d) / (n^2 * (n - 1))
            V <- V + f * (diag(at.risk * n) - outer(at.risk, at.risk))
        }
    }
    idx <- seq_len(k - 1)
    diff <- (O - E)[idx]
    chisq <- as.numeric(t(diff) %*% solve(V[idx, idx, drop = FALSE]) %*% diff)
    list(chisq = chisq, df = k - 1, p = stats::pchisq(chisq, k - 1, lower.tail = FALSE),
         observed = stats::setNames(O, lv), expected = stats::setNames(E, lv),
         n = stats::setNames(as.numeric(table(factor(group, levels = lv))), lv))
}

#' Split values at the median
#'
#' @param x Numeric values
#' @param labels Labels for below / at-or-above the median
#' @return Character vector (\code{NA} where \code{x} is \code{NA})
#' @export
median_split <- function(x, labels = c("low", "high")) {
    m <- stats::median(x, na.rm = TRUE)
    ifelse(is.na(x), NA_character_, ifelse(x >= m, labels[2], labels[1]))
}

#' Landmark survival status
#'
#' Survival status at a landmark time, e.g. alive at 1 year: subjects whose
#' time is at least \code{at} are \code{"alive"}; subjects with an event
#' before \code{at} are \code{"died"}; subjects censored before \code{at}
#' are \code{NA} (unknown).
#'
#' @param outcomes Data with time and event columns (e.g. \code{\link{subject_outcomes}})
#' @param time,event Column names
#' @param at Landmark time, in the units of \code{time} (months for PRINCE OS)
#' @param labels Labels for alive / died
#' @return Character vector, one per row of \code{outcomes}
#' @export
survival_status <- function(outcomes, time = "os", event = "os_event", at = 12,
                            labels = c(alive = "alive at 1 year", died = "died within 1 year")) {
    t <- outcomes[[time]]; e <- as.logical(outcomes[[event]])
    ifelse(is.na(t), NA_character_,
           ifelse(t >= at, labels[["alive"]], ifelse(!is.na(e) & e, labels[["died"]], NA_character_)))
}

#' Survival association of a subject-level biomarker
#'
#' Splits subjects at the median of \code{values} and tests the difference in
#' survival (log-rank). The workhorse behind "OS stratified by baseline
#' biomarker above/below the median" analyses.
#'
#' @param values Named numeric vector, names = subject ids (e.g. baseline
#'   frequency of a cell population)
#' @param outcomes Output of \code{\link{subject_outcomes}}
#' @param time,event Outcome columns
#' @return \code{outcomes} restricted to subjects with a value, plus
#'   \code{value} and \code{group} (\code{"low"}/\code{"high"}); the
#'   log-rank test is attached as attribute \code{"logrank"}
#' @export
survival_by_median <- function(values, outcomes, time = "os", event = "os_event") {
    values <- values[!is.na(values)]
    tab <- outcomes[outcomes$subject_id %in% names(values), , drop = FALSE]
    tab$value <- unname(values[tab$subject_id])
    tab$group <- median_split(tab$value)
    attr(tab, "logrank") <- logrank_test(tab[[time]], tab[[event]], tab$group)
    tab
}

#' Change from baseline per subject
#'
#' For each subject (and target, e.g. cell population or protein), the
#' change of \code{value} relative to the subject's value at the baseline
#' timepoint.
#'
#' @param tab Long data with \code{subject_id}, \code{timepoint_id},
#'   \code{value} (e.g. from \code{\link{add_sample_context}})
#' @param baseline Baseline timepoint id (e.g. \code{"C1D1"})
#' @param by Additional columns identifying a series (default: all target
#'   columns, i.e. columns other than the standard sample/context columns)
#' @param method \code{"log2_ratio"} (\code{log2(value / baseline)}; for
#'   frequencies), \code{"difference"} (\code{value - baseline}; for values
#'   already on a log scale, like Olink NPX) or \code{"ratio"}
#' @param pseudocount Added to both before a ratio
#' @return \code{tab} restricted to subjects with a baseline value, plus
#'   \code{baseline_value} and \code{change}
#' @export
change_from_baseline <- function(tab, baseline = "C1D1", by = NULL,
                                 method = c("log2_ratio", "difference", "ratio"), pseudocount = 0) {
    method <- match.arg(method)
    if (is.null(by)) by <- intersect(c("cell_population", "epitope_id", "hgnc_symbol", "signature",
                                       "measurement_set", "target"), names(tab))
    key <- do.call(paste, c(tab[c("subject_id", by)], sep = "\r"))
    base <- tab[tab$timepoint_id == baseline & !is.na(tab$value), , drop = FALSE]
    bkey <- do.call(paste, c(base[c("subject_id", by)], sep = "\r"))
    bval <- tapply(base$value, bkey, mean)
    tab$baseline_value <- unname(bval[key])
    tab <- tab[!is.na(tab$baseline_value), , drop = FALSE]
    tab$change <- switch(method,
                         log2_ratio = log2((tab$value + pseudocount) / (tab$baseline_value + pseudocount)),
                         difference = tab$value - tab$baseline_value,
                         ratio = (tab$value + pseudocount) / (tab$baseline_value + pseudocount))
    tab$change[!is.finite(tab$change)] <- NA
    tab
}
