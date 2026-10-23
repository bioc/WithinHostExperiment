# R/qc-temporal.R
# Longitudinal QC: flag transient iSNVs across timepoints

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames colData
#' @importFrom S4Vectors DataFrame
NULL

#' Flag temporally inconsistent iSNVs in longitudinal data
#'
#' For longitudinal samples from the same host, identifies iSNVs
#' that appear at only one timepoint (transient) versus those
#' detected at multiple consecutive timepoints (persistent).
#' Transient iSNVs are more likely to be sequencing artefacts;
#' persistent ones are more likely genuine within-host variants.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with longitudinal
#'   samples from the same host.
#' @param hostCol Character scalar. Column name in \code{colData}
#'   for host identifier (default: \code{"host_id"}).
#' @param timeCol Character scalar. Column name in \code{colData}
#'   for time point (default: \code{"timepoint"}). It must hold numbers,
#'   dates, date-times, time differences or a factor whose levels are in
#'   time order; text is read as numbers when every value is one, and is
#'   otherwise an error, since text such as \code{"day10"} and
#'   \code{"day2"} has no time order. A blank or \code{NA} is a missing
#'   timepoint.
#' @param minTimepoints Integer scalar. Minimum number of distinct
#'   timepoints at which an iSNV must be detected to be considered
#'   persistent (default: 2). Replicate samples taken at one
#'   timepoint count as one timepoint.
#' @param minConsecutive Integer scalar. Minimum number of
#'   \strong{consecutive} timepoints required, counted over the
#'   distinct timepoints at which the host was sampled. If \code{NULL}
#'   (default), any \code{minTimepoints} detections suffice
#'   regardless of consecutiveness. Set to 2 or 3 for stricter
#'   filtering.
#' @param usePassedOnly Logical. If \code{TRUE} (default), only
#'   consider QC-passed entries; a cell whose \code{qcPass} is
#'   \code{NA} counts as not passed.
#'
#' @return A \code{\link{WithinHostExperiment}} with updated
#'   \code{qcPass} assay and a new column \code{temporal_class} in
#'   \code{mcols(rowRanges)}. Sites are classified within each host;
#'   a transient site is flagged only in the samples of the hosts
#'   where it is transient. The site-level \code{temporal_class} is
#'   \code{"persistent"} if the site is persistent in at least one
#'   host, \code{"transient"} if it is transient in every host where
#'   it was detected, and \code{NA} where it was never detected.
#'
#' @details
#' This check applies only where a host was sampled at more than one
#' timepoint. Detections are counted per distinct timepoint: a site
#' called in any sample taken at a timepoint, including any of its
#' technical replicates, is detected at that timepoint, so replicates
#' alone cannot make a site persistent. A host sampled at a single
#' timepoint, however many replicates it has, is not classified.
#' Samples whose timepoint is \code{NA} are not counted, with a
#' warning, but are flagged with the rest of their host. Transient
#' sites are flagged by setting \code{qcPass} to \code{FALSE} where
#' it was \code{TRUE}; \code{NA} cells are left as they are.
#'
#' The biological rationale:
#' genuine within-host variants arise from viral replication and
#' should persist (or change smoothly) across consecutive
#' timepoints. An iSNV that appears once at 5\% and is absent at
#' all other timepoints is more likely a sequencing or library
#' preparation artefact.
#'
#' The function records the filtering step in \code{qcLog} with
#' step name \code{"temporal_consistency"}.
#'
#' @references
#' Farjo M et al. (2024). Within-host evolutionary dynamics and
#' tissue compartmentalization during acute SARS-CoV-2 infection.
#' \emph{J. Virol.} 98:e01618-23.
#' \doi{10.1128/jvi.01618-23}
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
#' freq <- matrix(c(0.1, NA,  0.2,   # var1: 2 timepoints
#'                  0.05, NA,  NA,    # var2: 1 timepoint (transient)
#'                  0.15, 0.2, 0.1),  # var3: 3 timepoints
#'                nrow = 3, ncol = 3)
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freq),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("t1", "t2", "t3"),
#'         host_id = c("H1", "H1", "H1"),
#'         timepoint = c(1, 2, 3)))
#' whe <- flagTemporalInconsistency(whe, minTimepoints = 2L)
#' S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
flagTemporalInconsistency <- function(whe,
                                      hostCol = "host_id",
                                      timeCol = "timepoint",
                                      minTimepoints = 2L,
                                      minConsecutive = NULL,
                                      usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    cd <- colData(whe)
    if (!hostCol %in% colnames(cd))
        stop("Column '", hostCol, "' not found in colData.")
    if (!timeCol %in% colnames(cd))
        stop("Column '", timeCol, "' not found in colData.")

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        matrix(TRUE, nrow = nrow(whe), ncol = ncol(whe))
    }

    hosts <- as.character(cd[[hostCol]])
    times <- .timeOrder(cd[[timeCol]], timeCol)
    n_sites <- nrow(freq_mat)

    no_time <- !is.na(hosts) & is.na(times)
    if (any(no_time))
        warning(sum(no_time), " sample(s) with NA in '", timeCol,
                "' are not counted towards persistence.")

    ## A cell is a detection when it has a frequency and passed QC; an
    ## NA qcPass cell counts as not passed
    det_mat <- !is.na(freq_mat) & !is.na(qc_mat) & qc_mat

    total_flagged <- 0L

    ## Group by host
    unique_hosts <- unique(hosts[!is.na(hosts)])
    class_by_host <- matrix(NA_character_, nrow = n_sites,
                            ncol = length(unique_hosts),
                            dimnames = list(NULL, unique_hosts))

    for (h in unique_hosts) {
        h_idx <- which(hosts == h & !is.na(times))
        h_times <- sort(unique(times[h_idx]))
        if (length(h_times) < 2L) next

        ## Collapse replicate samples to distinct timepoints, in time
        ## order: a site is detected at a timepoint if it is detected in
        ## any sample taken then
        det_t <- matrix(vapply(h_times, function(t) {
            cols <- h_idx[times[h_idx] == t]
            rowSums(det_mat[, cols, drop = FALSE]) > 0
        }, logical(n_sites)), nrow = n_sites)

        n_detected <- rowSums(det_t)
        is_persistent <- n_detected >= minTimepoints

        if (!is.null(minConsecutive)) {
            ## Longest run of consecutive timepoints with a detection
            run <- max_run <- integer(n_sites)
            for (k in seq_along(h_times)) {
                run <- ifelse(det_t[, k], run + 1L, 0L)
                max_run <- pmax(max_run, run)
            }
            is_persistent <- is_persistent & max_run >= minConsecutive
        }

        present <- n_detected > 0
        class_by_host[present, h] <- ifelse(is_persistent[present],
                                            "persistent", "transient")
    }

    ## Flag transient iSNVs in the samples of the host where they are
    ## transient; a site can be transient in one host and persistent in
    ## another. Cells with an NA qcPass keep it.
    qc <- assay(whe, "qcPass")
    for (h in unique_hosts) {
        transient_rows <- which(class_by_host[, h] == "transient")
        if (!length(transient_rows)) next
        host_cols <- which(hosts == h & !is.na(hosts))
        cells <- qc[transient_rows, host_cols, drop = FALSE]
        passed <- !is.na(cells) & cells
        cells[passed] <- FALSE
        qc[transient_rows, host_cols] <- cells
        total_flagged <- total_flagged + sum(passed)
    }
    assay(whe, "qcPass") <- qc

    ## Site-level summary across hosts
    temporal_class <- apply(class_by_host, 1L, function(classes) {
        classes <- classes[!is.na(classes)]
        if (!length(classes)) {
            NA_character_
        } else if (any(classes == "persistent")) {
            "persistent"
        } else {
            "transient"
        }
    })

    ## Add temporal_class to rowRanges
    rr <- SummarizedExperiment::rowRanges(whe)
    S4Vectors::mcols(rr)$temporal_class <- temporal_class
    SummarizedExperiment::rowRanges(whe) <- rr

    ## Log
    whe <- .addQcStep(whe, "temporal_consistency", "minTimepoints",
                      minTimepoints, total_flagged,
                      sum(qc, na.rm = TRUE))

    whe
}

## The timepoints of colData as numbers in time order: numbers, dates,
## date-times, time differences and factor level codes as they are, text
## only when every value is a number; a blank or an all-NA column is a
## missing timepoint
.timeOrder <- function(x, timeCol) {
    if (inherits(x, "POSIXlt")) x <- as.POSIXct(x)
    if (inherits(x, "difftime")) return(as.numeric(x, units = "secs"))
    if (is.factor(x) || is.numeric(x) || inherits(x, c("Date", "POSIXct")))
        return(as.numeric(x))
    ## NA alone is stored as logical
    if (is.logical(x) && all(is.na(x))) return(rep(NA_real_, length(x)))
    if (is.character(x)) {
        x[!is.na(x) & !nzchar(trimws(x))] <- NA
        number <- paste0("^\\s*[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)",
                         "([eE][-+]?[0-9]+)?\\s*$")
        bad <- !is.na(x) & !grepl(number, x)
        if (!any(bad)) return(as.numeric(x))
        stop("Column '", timeCol, "' holds text such as '", x[bad][1L],
             "', which has no time order. Give the timepoints as numbers, ",
             "dates or a factor with its levels in time order, for example ",
             "factor(x, levels = c(\"day2\", \"day10\")).", call. = FALSE)
    }
    stop("Column '", timeCol, "' must hold numbers, dates, date-times, ",
         "time differences or a factor, not ", class(x)[1L], ".",
         call. = FALSE)
}
