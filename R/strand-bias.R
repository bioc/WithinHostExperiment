# R/strand-bias.R
# Strand coverage imbalance (SCI) score for strand bias

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames
#' @importFrom S4Vectors DataFrame
NULL

# ============================================================
# sciStrandBias -- vectorised SCI computation
# ============================================================

#' Strand Coverage Imbalance (SCI) score
#'
#' Computes a strand bias score that normalises the per-site
#' forward/reverse alternative-allele ratio by the global strand
#' coverage ratio, making it suitable for amplicon data where strand
#' coverage is inherently unequal. Mostefai et al. (2024) showed that
#' strand bias filters which ignore this imbalance can discard genuine
#' variants; their own metric is a per-strand binomial likelihood.
#'
#' @param fwd_alt Numeric vector. Forward-strand alternative allele
#'   read counts. All four count vectors must be numeric and hold
#'   non-negative, finite values.
#' @param rev_alt Numeric vector. Reverse-strand alternative allele
#'   read counts.
#' @param fwd_total Numeric vector. Forward-strand total read counts.
#' @param rev_total Numeric vector. Reverse-strand total read counts.
#'
#' @return A numeric vector of SCI scores (log2 scale). Values
#'   near zero indicate no strand bias; large absolute values
#'   indicate strand-specific enrichment or depletion.
#'
#' @details
#' The SCI is defined as:
#' \deqn{\mathrm{SCI} = \log_2 \!\left(
#'   \frac{f_{\mathrm{alt}} / r_{\mathrm{alt}}}
#'        {F_{\mathrm{total}} / R_{\mathrm{total}}}
#' \right)}{SCI = log2((fwd_alt/rev_alt) / (sum(fwd_total)/sum(rev_total)))}
#'
#' where \eqn{F_{\mathrm{total}}} and \eqn{R_{\mathrm{total}}} are
#' the sums of forward and reverse total counts across all sites.
#' A pseudocount of 1 is added to strand counts to avoid division
#' by zero.
#'
#' @references
#' Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host
#' variation by leveraging large-scale sequencing data.
#' \emph{NAR Genomics and Bioinformatics} 6:lqae145.
#' \doi{10.1093/nargab/lqae145}
#'
#' @export
#' @examples
#' sci <- sciStrandBias(
#'     fwd_alt   = c(50, 2, 30),
#'     rev_alt   = c(45, 40, 28),
#'     fwd_total = c(500, 500, 500),
#'     rev_total = c(480, 480, 480)
#' )
#' sci
sciStrandBias <- function(fwd_alt, rev_alt, fwd_total, rev_total) {
    if (length(fwd_alt) != length(rev_alt) ||
        length(fwd_alt) != length(fwd_total) ||
        length(fwd_alt) != length(rev_total))
        stop("All input vectors must have the same length.")
    counts <- list(fwd_alt = fwd_alt, rev_alt = rev_alt,
                   fwd_total = fwd_total, rev_total = rev_total)
    for (nm in names(counts)) {
        x <- counts[[nm]]
        if (!is.numeric(x))
            stop("'", nm, "' must be numeric.")
        if (any(is.nan(x) | (!is.na(x) & (is.infinite(x) | x < 0))))
            stop("'", nm, "' must hold non-negative, finite counts.")
    }

    ## Pseudocount to avoid division by zero
    fwd_alt_s  <- fwd_alt + 1
    rev_alt_s  <- rev_alt + 1
    fwd_total_s <- fwd_total + 1
    rev_total_s <- rev_total + 1

    observed_ratio <- fwd_alt_s / rev_alt_s
    global_ratio   <- sum(fwd_total_s) / sum(rev_total_s)
    log2(observed_ratio / global_ratio)
}


# ============================================================
# flagStrandBiasSCI -- flag sites using SCI threshold
# ============================================================

#' Flag strand-biased iSNVs using the SCI metric
#'
#' Applies Strand Coverage Imbalance filtering to a
#' \code{\link{WithinHostExperiment}}. Sites whose absolute SCI
#' score exceeds \code{maxSCI} are flagged in the \code{qcPass}
#' assay.
#'
#' @param whe A \code{\link{WithinHostExperiment}} object. Must
#'   contain the assays \code{fwdAltCount}, \code{revAltCount},
#'   \code{fwdTotalCount}, and \code{revTotalCount}.
#' @param maxSCI Numeric scalar. Maximum allowable absolute SCI
#'   score (default: 1.5). Sites with \code{|SCI| > maxSCI} are
#'   flagged.
#'
#' @return A \code{\link{WithinHostExperiment}} with updated
#'   \code{qcPass} and \code{qcLog}.
#'
#' @details
#' Unlike the simple fold-difference strand bias in
#' \code{\link{flagISNV}}, this function uses the SCI score, which
#' normalises by the global strand coverage ratio (see
#' \code{\link{sciStrandBias}}). This is preferred for amplicon-based sequencing
#' where forward and reverse strand coverage is inherently
#' unbalanced.
#'
#' @references
#' Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host
#' variation by leveraging large-scale sequencing data.
#' \emph{NAR Genomics and Bioinformatics} 6:lqae145.
#' \doi{10.1093/nargab/lqae145}
#'
#' @seealso \code{\link{sciStrandBias}} for the underlying
#'   computation, \code{\link{flagISNV}} for standard QC.
#'
#' @export
#' @examples
#' library(SummarizedExperiment)
#' library(GenomicRanges)
#' gr <- GRanges("seg", IRanges::IRanges(start = 1:3, width = 1))
#' assays <- list(
#'     altFreq      = matrix(c(0.1, 0.05, 0.2), ncol = 1),
#'     totalDepth   = matrix(c(1000L, 1000L, 1000L), ncol = 1),
#'     qcPass       = matrix(TRUE, nrow = 3, ncol = 1),
#'     fwdAltCount  = matrix(c(50L, 2L, 30L), ncol = 1),
#'     revAltCount  = matrix(c(45L, 40L, 28L), ncol = 1),
#'     fwdTotalCount = matrix(c(500L, 500L, 500L), ncol = 1),
#'     revTotalCount = matrix(c(480L, 480L, 480L), ncol = 1)
#' )
#' whe <- WithinHostExperiment(
#'     rowRanges = gr, assays = assays,
#'     colData = S4Vectors::DataFrame(sample_id = "s1"))
#' whe <- flagStrandBiasSCI(whe, maxSCI = 1.5)
#' qcLog(whe)
flagStrandBiasSCI <- function(whe, maxSCI = 1.5) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    needed <- c("fwdAltCount", "revAltCount",
                 "fwdTotalCount", "revTotalCount")
    missing <- setdiff(needed, assayNames(whe))
    if (length(missing) > 0L)
        stop("Missing required assays: ",
             paste(missing, collapse = ", "))

    if (!"qcPass" %in% assayNames(whe))
        stop("WithinHostExperiment must have a 'qcPass' assay.")

    qc <- assay(whe, "qcPass")
    fwd_alt   <- assay(whe, "fwdAltCount")
    rev_alt   <- assay(whe, "revAltCount")
    fwd_total <- assay(whe, "fwdTotalCount")
    rev_total <- assay(whe, "revTotalCount")

    ## Compute SCI per column (sample)
    sci_mat <- matrix(NA_real_, nrow = nrow(qc), ncol = ncol(qc))
    for (j in seq_len(ncol(qc))) {
        sci_mat[, j] <- sciStrandBias(
            fwd_alt[, j], rev_alt[, j],
            fwd_total[, j], rev_total[, j]
        )
    }

    fail <- !is.na(sci_mat) & abs(sci_mat) > maxSCI & qc
    n_flag <- sum(fail, na.rm = TRUE)
    qc[fail] <- FALSE
    assay(whe, "qcPass") <- qc
    whe <- .addQcStep(whe, "strand_bias_sci_filter", "maxSCI",
                      maxSCI, n_flag,
                      sum(qc, na.rm = TRUE))
    whe
}
