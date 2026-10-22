# R/validation.R
# Validation and benchmarking utilities

#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
NULL

# ============================================================
# Heuristic plausibility check
# ============================================================

#' Check diversity estimates against heuristic plausibility ranges
#'
#' Compares observed diversity metrics against wide plausibility
#' ranges for common pathogens, as a check on units and orders of
#' magnitude. The ranges are heuristic: they are the package's own,
#' set generously enough that a correctly computed estimate on a
#' typical data set falls inside them, and they are not taken from
#' the literature or from any study. Each is listed in Details. A value
#' inside its range is therefore no evidence that an estimate is
#' correct, and one outside it most often points to a unit, depth or
#' QC problem rather than to biology.
#'
#' @details
#' The ranges, per site for \eqn{\pi} and on the scale
#' \code{\link{calcDiversity}} returns for Shannon entropy:
#' \tabular{lll}{
#'   \strong{pathogen} \tab \strong{pi} \tab \strong{shannon} \cr
#'   \code{sars_cov_2} \tab 1e-6 to 1e-3 \tab 0 to 15 \cr
#'   \code{influenza} \tab 1e-5 to 5e-3 \tab 0 to 20 \cr
#'   \code{hiv} \tab 1e-4 to 5e-2 \tab 0 to 30 \cr
#'   \code{tb} \tab 1e-9 to 1e-4 \tab 0 to 5 \cr
#'   \code{generic} \tab 0 to 1 \tab 0 to 100
#' }
#'
#' A wide table, as \code{\link{calcDiversity}} returns, is checked
#' row by row, one row per sample and metric; a long table with
#' \code{metric} and \code{value} columns is checked as given. Metrics
#' without a range are left out.
#'
#' @param div_df A \code{DataFrame} or \code{data.frame} with columns
#'   \code{metric} and \code{value}, or containing diversity columns
#'   such as \code{pi} and \code{shannon} (as returned by
#'   \code{\link{calcDiversity}}), with an optional \code{sample_id}.
#' @param pathogen Character scalar. Pathogen whose plausibility
#'   range is used.
#'   One of \code{"sars_cov_2"}, \code{"influenza"}, \code{"hiv"},
#'   \code{"tb"}, or \code{"generic"}.
#'
#' @return A \code{DataFrame} with one row per sample and metric and
#'   columns:
#'   \describe{
#'     \item{sample_id}{Sample identifier, taken from \code{div_df}
#'       when it has that column, otherwise the row number; \code{NA}
#'       for long input.}
#'     \item{metric}{Diversity metric name.}
#'     \item{observed}{Observed value.}
#'     \item{expected_min}{Lower bound of the plausibility range.}
#'     \item{expected_max}{Upper bound of the plausibility range.}
#'     \item{in_range}{Logical; whether the observed value falls
#'       within the range.}
#'   }
#'
#' @export
#' @examples
#' div <- S4Vectors::DataFrame(sample_id = c("a", "b"),
#'     pi = c(1e-4, 0.3), shannon = c(5.0, 7.5))
#' checkDiversityPlausibility(div, pathogen = "sars_cov_2")
checkDiversityPlausibility <- function(div_df,
                                       pathogen = c("sars_cov_2",
                                                    "influenza", "hiv",
                                                    "tb", "generic")) {
    pathogen <- match.arg(pathogen)

    ranges <- list(
        sars_cov_2 = list(pi = c(1e-6, 1e-3), shannon = c(0, 15)),
        influenza  = list(pi = c(1e-5, 5e-3), shannon = c(0, 20)),
        hiv        = list(pi = c(1e-4, 5e-2), shannon = c(0, 30)),
        tb         = list(pi = c(1e-9, 1e-4), shannon = c(0, 5)),
        generic    = list(pi = c(0, 1),       shannon = c(0, 100))
    )

    ref <- ranges[[pathogen]]

    ## Determine input format: long (metric/value cols) or wide
    if (is.data.frame(div_df) || is(div_df, "DataFrame")) {
        cn <- colnames(div_df)
        if (all(c("metric", "value") %in% cn)) {
            ## Long format: one row per metric, no sample
            metrics <- as.character(div_df$metric)
            values  <- as.numeric(div_df$value)
            samples <- rep(NA_character_, length(metrics))
        } else {
            ## Wide format, as calcDiversity() returns: every row is a
            ## sample, and every sample is checked
            avail <- intersect(cn, names(ref))
            if (length(avail) == 0L)
                stop("No recognised diversity metrics found in 'div_df'. ",
                     "Expected columns: ",
                     paste(names(ref), collapse = ", "))
            n <- nrow(div_df)
            ids <- if ("sample_id" %in% cn) as.character(div_df$sample_id)
                   else as.character(seq_len(n))
            metrics <- rep(avail, each = n)
            samples <- rep(ids, times = length(avail))
            values  <- unlist(lapply(avail, function(m)
                as.numeric(div_df[[m]])), use.names = FALSE)
        }
    } else {
        stop("'div_df' must be a DataFrame or data.frame.")
    }

    ## Filter to metrics with known ranges
    keep <- metrics %in% names(ref)
    metrics <- metrics[keep]
    values  <- values[keep]
    samples <- samples[keep]

    if (length(metrics) == 0L)
        stop("No metrics with known ranges for pathogen '", pathogen, "'.")

    exp_min <- vapply(metrics, function(m) ref[[m]][1L], numeric(1))
    exp_max <- vapply(metrics, function(m) ref[[m]][2L], numeric(1))
    in_range <- !is.na(values) & values >= exp_min & values <= exp_max

    DataFrame(
        sample_id    = samples,
        metric       = unname(metrics),
        observed     = unname(values),
        expected_min = unname(exp_min),
        expected_max = unname(exp_max),
        in_range     = unname(in_range)
    )
}

# ============================================================
# Caller benchmarking
# ============================================================

#' Benchmark variant callers against truth or each other
#'
#' If a truth \code{WithinHostExperiment} is provided, computes
#' true positives, false positives, false negatives, precision,
#' recall, and F1 for each caller. If no truth is provided,
#' computes pairwise agreement using
#' \code{\link{calcCallerConcordance}}.
#'
#' @param whe_list A named list of \code{\link{WithinHostExperiment}}
#'   objects, one per variant caller.
#' @param truth An optional \code{\link{WithinHostExperiment}}
#'   containing the true variant set. Default \code{NULL}.
#'
#' @return A \code{DataFrame}. When \code{truth} is provided:
#'   \describe{
#'     \item{caller}{Caller name.}
#'     \item{n_total}{Total sites called by this caller.}
#'     \item{TP}{True positives.}
#'     \item{FP}{False positives.}
#'     \item{FN}{False negatives.}
#'     \item{precision}{TP / (TP + FP).}
#'     \item{recall}{TP / (TP + FN).}
#'     \item{F1}{Harmonic mean of precision and recall.}
#'   }
#'   When \code{truth} is \code{NULL}, the result of
#'   \code{\link{calcCallerConcordance}} is returned.
#'
#' @export
#' @examples
#' gr1 <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr1)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr1)$alt <- c("T", "G", "A")
#' whe1 <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
#'     rowRanges = gr1,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' benchmarkCallers(list(caller1 = whe1), truth = whe1)
benchmarkCallers <- function(whe_list, truth = NULL) {
    if (!is.list(whe_list) || length(whe_list) < 1L)
        stop("'whe_list' must be a named list of ",
             "WithinHostExperiment objects.")
    if (is.null(names(whe_list)))
        stop("'whe_list' must be a named list (names = caller names).")
    are_whe <- vapply(whe_list, is, logical(1), "WithinHostExperiment")
    if (!all(are_whe))
        stop("All elements of 'whe_list' must be ",
             "WithinHostExperiment objects.")

    if (is.null(truth)) {
        if (length(whe_list) < 2L)
            stop("At least 2 callers required for pairwise agreement ",
                 "when no truth set is provided.")
        return(calcCallerConcordance(whe_list))
    }

    if (!is(truth, "WithinHostExperiment"))
        stop("'truth' must be a WithinHostExperiment object.")

    ## Build truth site keys
    truth_sites <- .siteKeys(truth)

    ## Evaluate each caller
    results <- lapply(names(whe_list), function(nm) {
        caller_sites <- .siteKeys(whe_list[[nm]])
        tp <- length(intersect(caller_sites, truth_sites))
        fp <- length(setdiff(caller_sites, truth_sites))
        fn <- length(setdiff(truth_sites, caller_sites))
        n_total <- length(caller_sites)
        prec   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
        rec    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
        f1     <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
            2 * prec * rec / (prec + rec)
        } else {
            NA_real_
        }
        data.frame(caller = nm, n_total = n_total,
                   TP = tp, FP = fp, FN = fn,
                   precision = prec, recall = rec, F1 = f1,
                   stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, results)
    DataFrame(out)
}

## Internal: extract chr:pos:alt site keys from a WHE
.siteKeys <- function(whe) {
    rr <- rowRanges(whe)
    mc <- mcols(rr)
    chr_vals <- as.character(seqnames(rr))
    pos_vals <- start(rr)
    alt_vals <- if ("alt" %in% colnames(mc)) {
        as.character(mc$alt)
    } else {
        rep(NA_character_, length(rr))
    }
    freq_mat <- assay(whe, "altFreq")
    has_data <- rowSums(!is.na(freq_mat)) > 0L
    paste(chr_vals[has_data], pos_vals[has_data],
          alt_vals[has_data], sep = ":")
}
