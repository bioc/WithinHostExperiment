# R/bottleneck-quick.R
# Lightweight beta-binomial bottleneck MLE

#' @include bridge.R
#' @importFrom stats dbinom pbeta dbeta optimize
NULL

#' Quick bottleneck size estimation (beta-binomial approximate MLE)
#'
#' Estimates transmission bottleneck size using the approximate
#' beta-binomial method of Sobel Leonard et al. (2017), which
#' does not model sequencing depth explicitly. For the exact
#' beta-binomial method (which accounts for finite read depth)
#' or other approaches (presence-absence, KL divergence,
#' Wright-Fisher), use the \pkg{ViralBottleneck} package
#' (Zheng et al. 2025, \emph{Virus Evolution} 11:veaf071)
#' via \code{\link{asViralBottleneckInput}}.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with transmission
#'   pairs defined in \code{colData}.
#' @param pairId Character. Which pair to estimate.
#' @param maxNb Integer. Maximum bottleneck size to evaluate
#'   (default: 200); a single whole number of at least 1.
#' @param threshold Numeric. Variant calling threshold (default: 0.03), a
#'   single number above 0 and below 0.5.
#'
#' @details
#' As in the published method, the donor variants used have frequencies
#' from \code{threshold} to \code{1 - threshold}; those above are left
#' out, with a message, since a donor variant at frequency 1 that the
#' recipient lacks has probability 0 under every \code{Nb}. The function
#' stops with an error when the log-likelihood is not finite at any
#' \code{Nb} of the grid, where it would otherwise return \code{Nb = 1}.
#'
#' @return A named list: Nb (MLE), ci (95\% CI), loglik,
#'   n_variants, n_excluded, pair_id. \code{n_excluded} is a named
#'   integer, as from \code{\link{exactBottleneck}}, of the donor
#'   variants at or above \code{threshold} left out:
#'   \code{recipient_failed_qc}, sites where QC removed a recipient call
#'   at or above \code{threshold}, which is not an absence (see
#'   \code{\link{asViralBottleneckInput}}); \code{unknown_depth}, always 0
#'   here, since this method does not use depths; and
#'   \code{donor_near_fixation}, donor variants above
#'   \code{1 - threshold}.
#'
#' @references
#' Sobel Leonard A et al. (2017). Transmission Bottleneck Size
#' Estimation from Pathogen Deep-Sequencing Data, with an Application to
#' Human Influenza A Virus. \emph{J Virol} 91:e00171-17.
#' \doi{10.1128/JVI.00171-17}
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' quickBottleneck(whe, pairId = "p1", maxNb = 50L)
quickBottleneck <- function(whe, pairId, maxNb = 200L, threshold = 0.03) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    maxNb <- .checkMaxNb(maxNb)
    threshold <- .checkThreshold(threshold)

    ## Get donor-recipient frequencies
    input <- .bottleneckInput(whe, pairId = pairId, threshold = threshold)
    kept <- .dropNearFixed(input$table, threshold, "quickBottleneck()")
    vb_input <- kept$table
    n_excluded <- .nExcluded(recipient_failed_qc = input$n_failed_qc,
                             donor_near_fixation = kept$n_near)
    if (nrow(vb_input) == 0L) {
        .emptyMessage(pairId, n_excluded)
        return(list(Nb = NA_real_, ci = c(NA_real_, NA_real_),
                    loglik = NA_real_, n_variants = 0L,
                    n_excluded = n_excluded, pair_id = pairId))
    }

    d_freqs <- vb_input$donor_freq
    r_freqs <- vb_input$recipient_freq

    ## Scan log-likelihood over Nb = 1 to maxNb
    ## Vectorized over variant sites for each Nb
    ll_vals <- vapply(seq_len(maxNb), function(nb) {
        site_lls <- vapply(seq_along(d_freqs), function(i) {
            .bbSiteLogLik(d_freqs[i], r_freqs[i], nb, threshold)
        }, numeric(1))
        sum(site_lls)
    }, numeric(1))
    .checkFiniteLogLik(ll_vals, maxNb, pairId)

    mle_idx <- which.max(ll_vals)
    nb_mle  <- mle_idx
    max_ll  <- ll_vals[mle_idx]

    ## 95% CI via likelihood ratio (chi-sq df=1, cutoff = 1.92)
    ci_cutoff <- max_ll - 1.92
    in_ci     <- which(ll_vals >= ci_cutoff)
    ci_lo     <- min(in_ci)
    ci_hi     <- max(in_ci)

    list(
        Nb         = as.numeric(nb_mle),
        ci         = c(as.numeric(ci_lo), as.numeric(ci_hi)),
        loglik     = max_ll,
        n_variants = nrow(vb_input),
        n_excluded = n_excluded,
        pair_id    = pairId
    )
}

## As in the published estimators that work from frequencies, use the
## donor variants from threshold to 1 - threshold: one at frequency 1 that
## the recipient lacks has probability 0 under every Nb
.dropNearFixed <- function(vb_input, threshold, what) {
    near_fixed <- vb_input$donor_freq > 1 - threshold
    n_near <- sum(near_fixed)
    if (n_near > 0L)
        message(what, " uses donor variants from threshold to ",
                "1 - threshold: ", n_near, " donor variant site(s) above ",
                "1 - threshold (", 1 - threshold, ") were left out.")
    list(table = vb_input[!near_fixed, , drop = FALSE], n_near = n_near)
}

## Stop when no Nb of the grid has a finite log-likelihood, where
## which.max() would return Nb = 1
.checkFiniteLogLik <- function(ll_vals, maxNb, pairId) {
    if (!any(is.finite(ll_vals)))
        stop("The log-likelihood is not finite at any Nb from 1 to ", maxNb,
             " for pair '", pairId, "', so there is no estimate: some ",
             "donor variant has probability 0 under every bottleneck ",
             "size, as one at frequency 1 that the recipient lacks has.",
             call. = FALSE)
    invisible(TRUE)
}

## Say that a pair has no donor variant to estimate from
.emptyMessage <- function(pairId, n_excluded) {
    if (sum(n_excluded) > 0L) {
        message("No donor variant above threshold is left for pair '",
                pairId, "' once those left out are set aside.")
    } else {
        message("No donor variants above threshold for pair '",
                pairId, "'.")
    }
}

#' @keywords internal
.bbSiteLogLik <- function(nu_d, nu_r, Nb, threshold, eps = 0.001) {
    present <- (nu_r >= threshold)
    log_probs <- numeric(Nb + 1L)

    for (k in 0:Nb) {
        log_p_k <- dbinom(k, Nb, nu_d, log = TRUE)
        if (k == 0L) {
            log_probs[k + 1L] <- if (!present) log_p_k else -Inf
        } else if (k == Nb) {
            log_probs[k + 1L] <- if (present) log_p_k else -Inf
        } else {
            alpha  <- k + eps
            beta_p <- Nb - k + eps
            log_p_obs <- if (present) {
                dbeta(nu_r, alpha, beta_p, log = TRUE)
            } else {
                pbeta(threshold, alpha, beta_p, log.p = TRUE)
            }
            log_probs[k + 1L] <- log_p_k + log_p_obs
        }
    }
    .logSumExp(log_probs)
}

## threshold as a single number above 0 and below 0.5, or an error: at 0
## every recipient frequency counts as present, and from 0.5 the window
## from threshold to 1 - threshold holds no donor variant
.checkThreshold <- function(threshold) {
    ok <- is.numeric(threshold) && length(threshold) == 1L &&
        is.finite(threshold) && threshold > 0 && threshold < 0.5
    if (!ok)
        stop("'threshold' must be a single number above 0 and below 0.5, ",
             "got ", paste(format(threshold), collapse = ", "), ".",
             call. = FALSE)
    threshold
}

## maxNb as a single whole number of at least 1, or an error
.checkMaxNb <- function(maxNb) {
    ok <- is.numeric(maxNb) && length(maxNb) == 1L && is.finite(maxNb) &&
        maxNb >= 1 && maxNb == round(maxNb) &&
        maxNb <= .Machine$integer.max
    if (!ok)
        stop("'maxNb' must be a positive integer, got ",
             paste(format(maxNb), collapse = ", "), ".")
    as.integer(maxNb)
}

#' @keywords internal
.logSumExp <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(-Inf)
    mx <- max(x)
    mx + log(sum(exp(x - mx)))
}
