# R/bottleneck-exact.R
# Exact beta-binomial, Wright-Fisher, and presence-absence bottleneck methods

#' @include bridge.R
#' @importFrom stats dbinom pbeta dbeta optimize quantile
NULL

#' Exact bottleneck size estimation
#'
#' Estimates transmission bottleneck size using one of three methods:
#' \describe{
#'   \item{\code{exact_bb}}{Beta-binomial likelihood over the number
#'     of transmitted variant virions, which models finite recipient
#'     read depth, in the spirit of Sobel Leonard et al. (2017).
#'     Requires \code{altCount} and \code{totalDepth} assays. Two
#'     deviations from the published method are deliberate: when all
#'     or none of the transmitted virions carry the variant, the
#'     recipient read counts are modelled with a fixed sequencing
#'     error probability of \code{threshold / 2}, and the recipient
#'     detection threshold is not integrated out. The likelihood needs
#'     the recipient's read depth at each donor variant, so it cannot
#'     use a site where that depth is unknown: such a site would add
#'     the same constant at every \code{Nb}. This includes every site
#'     that the recipient's VCF file does not list, so with objects read
#'     from VCF files a donor variant the recipient lacks usually enters
#'     \code{exact_bb} not at all, and loss in transmission is then
#'     scored only by the other estimators. Such sites are left out,
#'     with a message, and counted in \code{n_excluded}, not in
#'     \code{n_variants}.}
#'   \item{\code{wright_fisher}}{Approximate score in which one
#'     generation of binomial drift is combined with a beta kernel of
#'     concentration \code{max(Nb, 2)} around the drifted frequency.
#'     It mixes a density for detected variants with probabilities for
#'     lost variants, so it is a heuristic rather than a likelihood;
#'     prefer \code{exact_bb} or \code{presence_absence} for
#'     inference.}
#'   \item{\code{presence_absence}}{Presence-absence likelihood using
#'     \eqn{P(\mathrm{detect} | Nb, \nu) =
#'       1 - (1-\nu)^{Nb}}{P(detect|Nb,freq) = 1-(1-freq)^Nb}.}
#' }
#'
#' @param whe A \code{\link{WithinHostExperiment}} with transmission
#'   pairs defined in \code{colData}.
#' @param pairId Character. Which pair to estimate.
#' @param maxNb Integer. Maximum bottleneck size to evaluate
#'   (default: 500); a single whole number of at least 1.
#' @param method Character. One of \code{"exact_bb"},
#'   \code{"wright_fisher"}, or \code{"presence_absence"}.
#' @param nboot Integer. Number of bootstrap replicates for
#'   confidence intervals (default: 0, uses a likelihood-ratio
#'   interval). Bootstrapping resamples variant sites and uses the
#'   random number generator, so set a seed for reproducible
#'   intervals.
#' @param threshold Numeric. Variant calling threshold (default: 0.03).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{Nb}}{Numeric. Maximum likelihood estimate of
#'     bottleneck size over the grid \code{1:maxNb}. A warning is
#'     issued when the estimate reaches \code{maxNb}, because the
#'     value is then censored by the grid.}
#'   \item{\code{ci}}{Numeric vector of length 2. 95\% confidence
#'     interval (likelihood-ratio or bootstrap percentile), also
#'     bounded by \code{maxNb}.}
#'   \item{\code{loglik}}{Numeric. Log-likelihood at the MLE; for
#'     \code{"wright_fisher"}, the heuristic score at its maximiser.}
#'   \item{\code{method}}{Character. Method used. For
#'     \code{"wright_fisher"} the \code{Nb} and \code{ci} are the
#'     maximiser and interval of a heuristic score, and a message says
#'     so when the method is run.}
#'   \item{\code{n_variants}}{Integer. Number of donor variants used.}
#'   \item{\code{n_excluded}}{Named integer. Donor variants above
#'     \code{threshold} that were left out, by reason:
#'     \code{recipient_failed_qc}, sites where QC removed the
#'     recipient's call, which is not an absence (see
#'     \code{\link{asViralBottleneckInput}}), left out by every method;
#'     and \code{unknown_depth}, sites of unknown recipient depth, which
#'     only \code{exact_bb} leaves out.}
#'   \item{\code{pair_id}}{Character. Pair identifier.}
#' }
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
#' exactBottleneck(whe, pairId = "p1", maxNb = 50L,
#'     method = "presence_absence")
exactBottleneck <- function(whe, pairId, maxNb = 500L,
                            method = c("exact_bb", "wright_fisher",
                                       "presence_absence"),
                            nboot = 0L, threshold = 0.03) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    method <- match.arg(method)
    maxNb  <- .checkMaxNb(maxNb)
    nboot  <- as.integer(nboot)
    if (threshold < 0 || threshold > 1)
        stop("'threshold' must be between 0 and 1, got ", threshold, ".")

    if (method == "wright_fisher")
        message("method = \"wright_fisher\" is a heuristic score, not a ",
                "likelihood; its Nb is not an estimate on a par with ",
                "exact_bb or presence_absence.")
    switch(method,
        exact_bb          = .exactBBEstimate(whe, pairId, maxNb,
                                             nboot, threshold),
        wright_fisher     = .wrightFisherEstimate(whe, pairId, maxNb,
                                                  nboot, threshold),
        presence_absence  = .presenceAbsenceEstimate(whe, pairId, maxNb,
                                                     nboot, threshold)
    )
}

# ------------------------------------------------------------------
# Method 1: exact beta-binomial (finite read depth)
# ------------------------------------------------------------------

#' @keywords internal
.exactBBEstimate <- function(whe, pairId, maxNb, nboot, threshold) {
    pair_data <- .extractPairData(whe, pairId, threshold,
                                  need_counts = TRUE)
    n_excluded <- .nExcluded(recipient_failed_qc = pair_data$n_failed_qc,
                             unknown_depth = pair_data$n_unknown_depth)
    if (pair_data$n == 0L)
        return(.emptyResult(pairId, "exact_bb", n_excluded))

    d_freqs   <- pair_data$donor_freq
    r_alt     <- pair_data$recipient_alt
    r_depth   <- pair_data$recipient_depth
    n_var     <- pair_data$n

    ## Scan log-likelihood over Nb = 1..maxNb
    ll_vals <- vapply(seq_len(maxNb), function(nb) {
        site_lls <- vapply(seq_len(n_var), function(i) {
            .exactBBSiteLogLik(d_freqs[i], r_alt[i], r_depth[i],
                               nb, threshold)
        }, numeric(1))
        sum(site_lls)
    }, numeric(1))

    .finalizeMLE(ll_vals, maxNb, nboot, n_var, pairId, "exact_bb",
                 n_excluded, d_freqs, r_alt, r_depth, threshold)
}

#' @keywords internal
.exactBBSiteLogLik <- function(nu_d, k_r, n_r, Nb, threshold,
                                eps = 0.001) {
    log_probs <- numeric(Nb + 1L)
    for (j in 0:Nb) {
        log_p_j <- dbinom(j, Nb, nu_d, log = TRUE)
        if (j == 0L) {
            ## All bottleneck virions are reference: recipient gets
            ## zero alt reads only if below threshold detection
            log_p_obs <- dbinom(k_r, n_r, threshold / 2, log = TRUE)
            log_probs[j + 1L] <- log_p_j + log_p_obs
        } else if (j == Nb) {
            ## All bottleneck virions are alt
            log_p_obs <- dbinom(k_r, n_r, 1 - threshold / 2,
                                log = TRUE)
            log_probs[j + 1L] <- log_p_j + log_p_obs
        } else {
            alpha  <- j + eps
            beta_p <- Nb - j + eps
            ## Expected recipient frequency drawn from Beta, then
            ## read counts from Binomial.  Integrate analytically
            ## via beta-binomial PMF.
            log_p_obs <- .logBetaBinomPMF(k_r, n_r, alpha, beta_p)
            log_probs[j + 1L] <- log_p_j + log_p_obs
        }
    }
    .logSumExp(log_probs)
}

#' Log beta-binomial PMF: log BetaBinom(k | n, alpha, beta)
#' @keywords internal
#' @noRd
.logBetaBinomPMF <- function(k, n, alpha, beta_p) {
    lchoose(n, k) +
        lbeta(k + alpha, n - k + beta_p) -
        lbeta(alpha, beta_p)
}

# ------------------------------------------------------------------
# Method 2: Wright-Fisher drift simulation
# ------------------------------------------------------------------

#' @keywords internal
.wrightFisherEstimate <- function(whe, pairId, maxNb, nboot,
                                  threshold) {
    input <- .bottleneckInput(whe, pairId = pairId, threshold = threshold)
    vb_input <- input$table
    n_excluded <- .nExcluded(recipient_failed_qc = input$n_failed_qc)
    if (nrow(vb_input) == 0L)
        return(.emptyResult(pairId, "wright_fisher", n_excluded))

    d_freqs <- vb_input$donor_freq
    r_freqs <- vb_input$recipient_freq
    present <- r_freqs >= threshold
    n_var   <- length(d_freqs)

    ## For each Nb, compute log-likelihood from one generation
    ## of binomial drift: k ~ Binom(Nb, nu_d), nu_r = k/Nb
    ll_vals <- vapply(seq_len(maxNb), function(nb) {
        site_lls <- vapply(seq_len(n_var), function(i) {
            .wfSiteLogLik(d_freqs[i], r_freqs[i], present[i],
                          nb, threshold)
        }, numeric(1))
        sum(site_lls)
    }, numeric(1))

    .finalizeMLE(ll_vals, maxNb, nboot, n_var, pairId,
                 "wright_fisher", n_excluded, d_freqs, r_freqs, present,
                 threshold)
}

#' @keywords internal
.wfSiteLogLik <- function(nu_d, nu_r, present, Nb, threshold,
                           eps = 1e-10) {
    log_probs <- numeric(Nb + 1L)
    for (k in 0:Nb) {
        log_p_k <- dbinom(k, Nb, nu_d, log = TRUE)
        drift_freq <- k / Nb
        if (!present) {
            ## Variant lost: drift frequency below threshold
            log_p_obs <- if (drift_freq < threshold) {
                0
            } else {
                -Inf
            }
        } else {
            ## Variant present: use density near the observed freq
            ## Approximate with narrow beta centred at drift_freq
            if (k == 0L) {
                log_p_obs <- -Inf
            } else if (k == Nb) {
                log_p_obs <- if (nu_r >= 1 - threshold) 0 else -Inf
            } else {
                ## Concentration proportional to Nb
                conc <- max(Nb, 2)
                a <- drift_freq * conc + eps
                b <- (1 - drift_freq) * conc + eps
                log_p_obs <- dbeta(nu_r, a, b, log = TRUE)
            }
        }
        log_probs[k + 1L] <- log_p_k + log_p_obs
    }
    .logSumExp(log_probs)
}

# ------------------------------------------------------------------
# Method 3: presence-absence
# ------------------------------------------------------------------

#' @keywords internal
.presenceAbsenceEstimate <- function(whe, pairId, maxNb, nboot,
                                     threshold) {
    input <- .bottleneckInput(whe, pairId = pairId, threshold = threshold)
    vb_input <- input$table
    n_excluded <- .nExcluded(recipient_failed_qc = input$n_failed_qc)
    if (nrow(vb_input) == 0L)
        return(.emptyResult(pairId, "presence_absence", n_excluded))

    d_freqs <- vb_input$donor_freq
    r_freqs <- vb_input$recipient_freq
    present <- r_freqs >= threshold
    n_var   <- length(d_freqs)

    ## Log-likelihood: product over sites of
    ##   present:  log(1 - (1 - nu_d)^Nb)
    ##   absent:   Nb * log(1 - nu_d)
    ll_vals <- vapply(seq_len(maxNb), function(nb) {
        site_lls <- vapply(seq_len(n_var), function(i) {
            .paSiteLogLik(d_freqs[i], present[i], nb)
        }, numeric(1))
        sum(site_lls)
    }, numeric(1))

    .finalizeMLE(ll_vals, maxNb, nboot, n_var, pairId,
                 "presence_absence", n_excluded, d_freqs, r_freqs, present,
                 threshold)
}

#' @keywords internal
.paSiteLogLik <- function(nu_d, present, Nb) {
    log_q <- Nb * log(1 - nu_d)
    if (present) {
        log(1 - exp(log_q))
    } else {
        log_q
    }
}

# ------------------------------------------------------------------
# Shared helpers
# ------------------------------------------------------------------

#' Extract donor/recipient data for a pair
#' @keywords internal
#' @noRd
.extractPairData <- function(whe, pairId, threshold,
                             need_counts = FALSE) {
    input <- .bottleneckInput(whe, pairId = pairId, threshold = threshold)
    vb_input <- input$table
    if (nrow(vb_input) == 0L)
        return(list(n = 0L, n_failed_qc = input$n_failed_qc,
                    n_unknown_depth = 0L))

    d_freqs <- vb_input$donor_freq
    n_var   <- length(d_freqs)

    if (!need_counts) {
        r_freqs <- vb_input$recipient_freq
        return(list(donor_freq = d_freqs, recipient_freq = r_freqs,
                    n = n_var, n_failed_qc = input$n_failed_qc))
    }

    ## Need per-site read counts for exact method
    anames <- assayNames(whe)
    if (!("altCount" %in% anames) || !("totalDepth" %in% anames))
        stop("Method 'exact_bb' requires 'altCount' and 'totalDepth' ",
             "assays in the WithinHostExperiment object.")

    ## The recipient's counts at the rows the donor variants come from,
    ## so that sites are matched by row, not by position, which is
    ## ambiguous for segmented genomes and multi-allelic sites
    idx <- cbind(input$sites$row_index, input$sites$recipient_col)
    r_alt   <- as.integer(as.matrix(assay(whe, "altCount"))[idx])
    r_depth <- as.integer(as.matrix(assay(whe, "totalDepth"))[idx])
    d_freqs <- as.numeric(d_freqs)

    ## A site of unknown recipient depth adds the same constant to the
    ## log-likelihood at every Nb, so it is left out and reported
    known <- !is.na(r_depth)
    n_unknown <- sum(!known)
    if (n_unknown > 0L)
        message("exact_bb uses no donor variant whose recipient depth is ",
                "unknown: ", n_unknown, " donor variant site(s) whose ",
                "recipient depth is not recorded were left out and are ",
                "counted in n_excluded, not in n_variants.")
    r_alt   <- r_alt[known]
    r_depth <- r_depth[known]
    d_freqs <- d_freqs[known]

    ## a missing alternative count at a known depth is taken as 0
    r_alt[is.na(r_alt)] <- 0L
    ## Ensure alt does not exceed depth
    r_alt <- pmin(r_alt, r_depth)

    list(donor_freq      = d_freqs,
         recipient_alt   = r_alt,
         recipient_depth = r_depth,
         n               = length(d_freqs),
         n_failed_qc     = input$n_failed_qc,
         n_unknown_depth = n_unknown)
}

#' Finalize MLE from log-likelihood vector
#' @keywords internal
#' @noRd
.finalizeMLE <- function(ll_vals, maxNb, nboot, n_var, pairId,
                          method, n_excluded, ...) {
    mle_idx <- which.max(ll_vals)
    nb_mle  <- mle_idx
    max_ll  <- ll_vals[mle_idx]
    if (nb_mle == maxNb)
        warning("The maximum likelihood estimate reached maxNb (", maxNb,
                "); the bottleneck size is censored by the grid. ",
                "Increase maxNb.", call. = FALSE)

    if (nboot > 0L) {
        ci <- .bootstrapCI(nboot, maxNb, n_var, method, ...)
    } else {
        ## 95% CI via likelihood ratio (chi-sq df=1, cutoff 1.92)
        ci_cutoff <- max_ll - 1.92
        in_ci     <- which(ll_vals >= ci_cutoff)
        ci <- c(min(in_ci), max(in_ci))
    }

    list(
        Nb         = as.numeric(nb_mle),
        ci         = as.numeric(ci),
        loglik     = max_ll,
        method     = method,
        n_variants = as.integer(n_var),
        n_excluded = n_excluded,
        pair_id    = pairId
    )
}

## Donor variants above threshold left out of an estimate, by reason
.nExcluded <- function(recipient_failed_qc = 0L, unknown_depth = 0L) {
    c(recipient_failed_qc = as.integer(recipient_failed_qc),
      unknown_depth = as.integer(unknown_depth))
}

#' Bootstrap confidence interval by resampling variant sites
#' @keywords internal
#' @noRd
.bootstrapCI <- function(nboot, maxNb, n_var, method, ...) {
    args <- list(...)
    boot_nbs <- vapply(seq_len(nboot), function(b) {
        idx <- sample.int(n_var, replace = TRUE)

        boot_ll <- vapply(seq_len(maxNb), function(nb) {
            site_lls <- vapply(idx, function(i) {
                if (method == "exact_bb") {
                    .exactBBSiteLogLik(args[[1]][i], args[[2]][i],
                                       args[[3]][i], nb, args[[4]])
                } else if (method == "wright_fisher") {
                    .wfSiteLogLik(args[[1]][i], args[[2]][i],
                                  args[[3]][i], nb, args[[4]])
                } else {
                    .paSiteLogLik(args[[1]][i], args[[3]][i], nb)
                }
            }, numeric(1))
            sum(site_lls)
        }, numeric(1))

        which.max(boot_ll)
    }, numeric(1))

    stats::quantile(boot_nbs, probs = c(0.025, 0.975),
                    names = FALSE)
}

#' Empty result when no variants found
#' @keywords internal
#' @noRd
.emptyResult <- function(pairId, method, n_excluded) {
    .emptyMessage(pairId, n_excluded)
    list(
        Nb         = NA_real_,
        ci         = c(NA_real_, NA_real_),
        loglik     = NA_real_,
        method     = method,
        n_variants = 0L,
        n_excluded = n_excluded,
        pair_id    = pairId
    )
}
