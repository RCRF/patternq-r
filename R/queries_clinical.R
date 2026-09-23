# Timepoints, clinical observations, adverse events, clinical interventions,
# and the cell population / TCR / microbiome entities of measurement sets.

#' Timepoints
#' @inheritParams samples
#' @return \code{timepoint_id}, \code{timepoint_relative_order},
#'   \code{timepoint_type}, \code{timepoint_offset}, cycle/day where present,
#'   ordered by relative order
#' @export
timepoints <- function(db = NULL, ...) {
    df <- do_query(dq(find = list(list("pull", "?t", list("*", list(":timepoint/type" = ":db/ident")))),
                      where = list(c("?t", ":timepoint/id"))), db = db, ...)
    if ("timepoint_relative_order" %in% names(df)) {
        df <- keep_provenance(df[order(df$timepoint_relative_order), , drop = FALSE], df)
        rownames(df) <- NULL
    }
    keep_provenance(order_columns(df, c("timepoint_id", "timepoint_relative_order", "timepoint_type")), df)
}

#' Clinical observation sets
#' @inheritParams samples
#' @return \code{clinical_observation_set_name}, \code{..._description}
#' @export
clinical_observation_sets <- function(db = NULL, ...) {
    do_query(dq(find = list(list("pull", "?c", c(":clinical-observation-set/name",
                                                 ":clinical-observation-set/description"))),
                where = list(c("?c", ":clinical-observation-set/name"))), db = db, ...)
}

#' Clinical observations
#'
#' Either all attributes of the observations in a set (one row per
#' observation), or one observation type across the dataset (one row per
#' subject/timepoint value).
#'
#' @inheritParams samples
#' @param obs.type Observation attribute without namespace, e.g. \code{"os"},
#'   \code{"pfs"}, \code{"bor"}, \code{"recist"}, \code{"ldh"}
#' @param set.name Clinical observation set name (alternative to
#'   \code{obs.type})
#' @param subjects Optional subject ids
#' @return With \code{obs.type}: \code{subject_id}, \code{timepoint_id},
#'   and a column named after the type (enum values as names). With
#'   \code{set.name}: one row per observation with all its attributes.
#' @export
clinical_observations <- function(db = NULL, obs.type = NULL, set.name = NULL, subjects = NULL, ...) {
    if (is.null(obs.type) && is.null(set.name))
        stop("Give obs.type or set.name", call. = FALSE)
    if (!is.null(obs.type)) {
        attr <- paste0(":clinical-observation/", sub("^:?clinical-observation/", "", obs.type))
        q <- dq(find = list("?subject-id", list("pull", "?o", list(list(":clinical-observation/timepoint" = ":timepoint/id"),
                                                                   list(":clinical-observation/study-day" = ":study-day/id"),
                                                                   attr))),
                where = list(c("?o", attr), c("?o", ":clinical-observation/subject", "?p"),
                             c("?p", ":subject/id", "?subject-id")),
                `in` = if (!is.null(subjects)) list(c("?subject-id", "...")),
                args = if (!is.null(subjects)) list(I(subjects)) else list())
        df <- do_query(q, db = db, ...)
        col <- clean_names(attr)
        names(df)[names(df) == col] <- clean_names(sub("^:?clinical-observation/", "", obs.type))
        return(df)
    }
    q <- dq(find = list(list("pull", "?o", list(
        "*",
        list(":clinical-observation/subject" = ":subject/id"),
        list(":clinical-observation/timepoint" = ":timepoint/id"),
        list(":clinical-observation/study-day" = ":study-day/id"),
        list(":clinical-observation/recist" = ":db/ident"),
        list(":clinical-observation/bor" = ":db/ident"),
        list(":clinical-observation/pfs-reason" = ":db/ident"),
        list(":clinical-observation/os-reason" = ":db/ident"),
        list(":clinical-observation/disease-stage" = ":db/ident"),
        list(":clinical-observation/metastasis-gdc-anatomic-sites" = ":gdc-anatomic-site/name")))),
        `in` = list("?set-name"), args = list(set.name),
        where = list(c("?cos", ":clinical-observation-set/name", "?set-name"),
                     c("?cos", ":clinical-observation-set/clinical-observations", "?o")))
    df <- do_query(q, db = db, ...)
    names(df) <- sub("^clinical_observation_", "", names(df))
    df
}

#' Subject outcomes: best overall response, PFS and OS
#'
#' One row per subject with \code{bor} (taken from \code{:clinical-observation/bor},
#' or derived from RECIST observations when absent: CR > PR > SD > PD),
#' \code{pfs}, \code{pfs_event}, \code{os}, \code{os_event} where present.
#' Errors if a subject has more than one value of a single-valued outcome.
#'
#' @inheritParams samples
#' @return \code{data.frame}
#' @export
subject_outcomes <- function(db = NULL, ...) {
    db <- ensure_db(db)
    get1 <- function(type) {
        r <- tryCatch(clinical_observations(db, obs.type = type, ...), error = function(e) NULL)
        if (is.null(r) || !nrow(r)) return(NULL)
        col <- clean_names(type)
        r <- r[!is.na(r[[col]]), c("subject_id", col), drop = FALSE]
        if (anyDuplicated(r$subject_id))
            stop(sprintf("More than one %s value for some subjects", type), call. = FALSE)
        r
    }
    bor <- get1("bor")
    if (is.null(bor)) {
        recist <- tryCatch(clinical_observations(db, obs.type = "recist", ...), error = function(e) NULL)
        if (!is.null(recist) && nrow(recist)) {
            rank <- c(CR = 1, PR = 2, SD = 3, PD = 4)
            best <- tapply(recist$recist, recist$subject_id, function(x) {
                x <- x[x %in% names(rank)]
                if (!length(x)) "Unknown" else names(rank)[min(rank[x])]
            })
            bor <- data.frame(subject_id = names(best), bor = unname(best), stringsAsFactors = FALSE)
        }
    }
    parts <- Filter(Negate(is.null), list(bor, get1("pfs"), get1("pfs-event"), get1("os"), get1("os-event")))
    ids <- data.frame(subject_id = subjects(db)$subject_id, stringsAsFactors = FALSE)
    out <- Reduce(function(a, b) merge(a, b, by = "subject_id", all.x = TRUE), parts, ids)
    with_provenance(out, db, provenance(ids)$basis_t)
}

#' Adverse events
#' @inheritParams samples
#' @param set.name Optional clinical observation set name
#' @return One row per adverse event
#' @export
adverse_events <- function(db = NULL, set.name = NULL, ...) {
    where <- list(c("?o", ":adverse-event/subject"))
    ins <- NULL; args <- list()
    if (!is.null(set.name)) {
        where <- list(c("?cos", ":clinical-observation-set/name", "?set-name"),
                      c("?cos", ":clinical-observation-set/adverse-events", "?o"))
        ins <- list("?set-name"); args <- list(set.name)
    }
    df <- do_query(dq(find = list(list("pull", "?o", list(
        "*",
        list(":adverse-event/subject" = ":subject/id"),
        list(":adverse-event/timepoint" = ":timepoint/id"),
        list(":adverse-event/meddra-adverse-event" = ":meddra-disease/preferred-name"),
        list(":adverse-event/ctcae-grade" = ":db/ident"),
        list(":adverse-event/ae-causality" = ":db/ident"),
        list(":adverse-event/study-day" = ":study-day/id")))),
        where = where, `in` = ins, args = args), db = db, ...)
    df
}

#' Clinical interventions (treatments, surgeries, biopsies, ...)
#' @inheritParams samples
#' @param subjects Optional subject ids
#' @return One row per intervention, with treatment regimen and drug names
#'   where present
#' @export
clinical_interventions <- function(db = NULL, subjects = NULL, ...) {
    q <- dq(find = list(list("pull", "?ci", list(
        "*",
        list(":clinical-intervention/subject" = ":subject/id"),
        list(":clinical-intervention/timepoint" = c(":timepoint/id", ":timepoint/relative-order")),
        list(":clinical-intervention/treatment-regimen" = list(
            ":treatment-regimen/name",
            list(":treatment-regimen/drug-regimens" = list(list(":drug-regimen/drug" = ":drug/preferred-name"),
                                                           ":drug-regimen/freetext-drug")))),
        list(":clinical-intervention/surgery-type" = ":db/ident"),
        list(":clinical-intervention/cancer-medication-category" = ":db/ident"),
        list(":clinical-intervention/radiation-therapy-category" = ":db/ident"),
        list(":clinical-intervention/biospecimen-type" = ":db/ident"),
        list(":clinical-intervention/biospecimen-collection" = ":db/ident"),
        list(":clinical-intervention/biospecimen-derived-samples" = ":sample/id")))),
        where = list(c("?ci", ":clinical-intervention/subject", "?p"), c("?p", ":subject/id", "?subject-id")),
        `in` = if (!is.null(subjects)) list(c("?subject-id", "...")),
        args = if (!is.null(subjects)) list(I(subjects)) else list())
    df <- do_query(q, db = db, ...)
    names(df) <- sub("^clinical_intervention_", "", names(df))
    df
}

ms_entities <- function(ref.attr, pattern, measurement.set, db, ...) {
    q <- dq(find = list(list("pull", "?e", pattern)),
            where = list(c("?ms", ":measurement-set/name", "?ms-name"), c("?ms", ref.attr, "?e")),
            `in` = list("?ms-name"), args = list(measurement.set))
    do_query(q, db = db, ...)
}

#' Entities of a measurement set
#'
#' Cell populations, TCRs, OTUs and SGBs (metagenomic species) attached to a
#' measurement set.
#'
#' @inheritParams samples
#' @param measurement.set Measurement set name
#' @return One row per entity
#' @name measurement_set_entities
NULL

#' @rdname measurement_set_entities
#' @export
cell_populations <- function(measurement.set, db = NULL, ...) {
    ms_entities(":measurement-set/cell-populations",
                list("*", list(":cell-population/cell-type" = ":cell-type/co-name"),
                     list(":cell-population/positive-markers" = ":epitope/id"),
                     list(":cell-population/negative-markers" = ":epitope/id"),
                     list(":cell-population/parent" = ":cell-population/name")),
                measurement.set, db, ...)
}

#' @rdname measurement_set_entities
#' @export
tcrs <- function(measurement.set, db = NULL, ...) {
    ms_entities(":measurement-set/tcrs", list("*"), measurement.set, db, ...)
}

#' @rdname measurement_set_entities
#' @export
otus <- function(measurement.set, db = NULL, ...) {
    ms_entities(":measurement-set/otus", list("*"), measurement.set, db, ...)
}

#' @rdname measurement_set_entities
#' @export
sgbs <- function(measurement.set, db = NULL, ...) {
    ms_entities(":measurement-set/sgbs", list("*"), measurement.set, db, ...)
}
