# Optional transit response formats for direct (uncached) queries: the query
# service answers POST /query with transit when asked (Accept
# application/transit+json or application/transit+msgpack). Decoding needs the
# optional transit package (remotes::install_github("vendekagon-labs/transit-r")).
# Decoded values are converted to what jsonlite gives for the JSON response,
# so everything downstream of raw_query works unchanged.

response_formats <- c(json = "application/json",
                      "transit+json" = "application/transit+json",
                      "transit+msgpack" = "application/transit+msgpack")

check_format <- function(format) {
    if (!is.character(format) || length(format) != 1 || !format %in% names(response_formats))
        stop("format must be one of: ", paste(names(response_formats), collapse = ", "), call. = FALSE)
    format
}

# leaf conversions (rapply walks the value in C, calling these only on matches)
transit_leaf_classes <- c("transit_keyword", "transit_symbol", "transit_int64", "transit_bigint",
                          "transit_decimal", "POSIXct", "transit_uuid", "transit_uri")

transit_leaf <- function(x) {
    if (inherits(x, c("transit_keyword", "transit_symbol"))) return(paste0(":", unclass(x)))
    if (inherits(x, c("transit_int64", "transit_bigint", "transit_decimal"))) return(as.numeric(unclass(x)))
    if (inherits(x, "POSIXct")) {
        # as the JSON response writes instants: UTC, milliseconds without trailing zeros
        # (from integer milliseconds: %OS3 truncates, e.g. .903 -> .902)
        ms <- round(as.numeric(x) * 1000)
        s <- format(as.POSIXct(ms %/% 1000, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%dT%H:%M:%S", tz = "UTC")
        frac <- sub("0+$", "", sprintf("%03d", as.integer(ms %% 1000)))
        return(paste0(s, ifelse(nzchar(frac), paste0(".", frac), ""), "Z"))
    }
    as.character(unclass(x))
}

# maps -> named lists, sets / lists / tagged values -> plain lists
transit_unmap <- function(x) {
    if (!is.list(x)) return(x)
    if (inherits(x, "transit_map")) {
        m <- unclass(x)
        out <- lapply(m$values, transit_unmap)
        names(out) <- vapply(m$keys, as.character, "")
        return(out)
    }
    if (is.object(x)) x <- unclass(x)
    nested <- vapply(x, is.list, NA)
    if (any(nested)) x[nested] <- lapply(x[nested], transit_unmap)
    x
}

#' Decode a transit query response into the JSON response's shape
#'
#' Used by \code{\link{raw_query}} for the \code{"transit+json"} and
#' \code{"transit+msgpack"} formats.
#'
#' @param body Raw response body
#' @param format \code{"transit+json"} or \code{"transit+msgpack"}
#' @return A list shaped like \code{jsonlite::fromJSON(json, simplifyVector = FALSE)}
#'   of the JSON response
#' @keywords internal
decode_transit <- function(body, format) {
    if (!requireNamespace("transit", quietly = TRUE))
        stop("format = \"", format, "\" needs the optional transit package: ",
             "remotes::install_github(\"vendekagon-labs/transit-r\")", call. = FALSE)
    x <- transit::from_transit(body, if (format == "transit+msgpack") "msgpack" else "json")
    x <- rapply(list(x), transit_leaf, classes = transit_leaf_classes, how = "replace")[[1]]
    transit_unmap(x)
}
