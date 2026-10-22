# R/neutrality.R
# Neutrality tests for within-host pathogen populations

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors DataFrame
#' @importFrom stats median pnorm
#' @importFrom methods setGeneric setMethod is
NULL

# ============================================================
# Generic
# ============================================================

#' Compute neutrality test statistics
#'
#' Generic for computing neutrality tests from within-host
#' variant data. Methods should return a \code{DataFrame} with
#' one row per sample.
#'
#' @param x A within-host variant container.
#' @param ... Arguments for methods.
#' @return A \code{DataFrame} of neutrality statistics.
#'
#' @export
#' @rdname calcNeutralityTests
setGeneric("calcNeutralityTests", function(x, ...)
    standardGeneric("calcNeutralityTests"))

# ============================================================
# Standalone math functions (exported, work on plain vectors)
# ============================================================

#' Tajima's D statistic
#'
#' Computes Tajima's D (Tajima 1989) from summary statistics of
#' a within-host pathogen population. The statistic measures the
#' difference between nucleotide diversity (\eqn{\hat{\pi}}) and
#' Watterson's estimator (\eqn{\hat{\theta}_W = S / a_1}),
#' normalised by its expected standard deviation under neutrality:
#' \deqn{D = \frac{\hat{\pi} - S/a_1}{\sqrt{e_1 S + e_2 S(S-1)}}}
#' A significantly negative D indicates an excess of rare
#' variants (e.g., population expansion or purifying selection),
#' while a positive D suggests balancing selection or population
#' contraction.
#'
#' @param S Integer scalar. Number of segregating sites.
#' @param n Numeric scalar. Sample size. For deep sequencing data,
#'   do \strong{not} use raw read depth directly -- the classical
#'   Tajima's D assumes independent lineages, whereas sequencing
#'   reads oversample a much smaller viral population.
#'   Use an effective sample size (e.g., \code{n = 50} or
#'   \code{n = 100}) that reflects the number of independently
#'   sampled viral genomes, or cap the depth at a biologically
#'   meaningful value. See Details.
#' @param pi_hat Numeric scalar. Nucleotide diversity estimate
#'   (per-site).
#' @param genomeLength Integer scalar. Reference genome length
#'   in bp.
#'
#' @return Numeric scalar. Tajima's D, or \code{NA} with a warning when
#'   \code{S} is below three or \code{n} rounds to fewer than four.
#'
#' @details
#' \strong{Fewer than three segregating sites:} with one segregating
#' site D is a monotone function of that site's frequency, and with two
#' it is a function of two frequencies; neither describes a frequency
#' spectrum, and at \code{S = 0} the statistic is 0/0. The function
#' therefore returns \code{NA}, with a warning, whenever \code{S < 3}.
#' Report \code{S} alongside every D, as
#' \code{\link{calcNeutralityTests}} does in its \code{n_sites}
#' column. With two or three sequences the variance of D is zero
#' (\eqn{c_1 = c_2 = 0}), so \code{n} must round to at least four;
#' below that D is likewise \code{NA}.
#'
#' \strong{Deep sequencing caveat:} The classical Tajima (1989)
#' formula was designed for independent haplotype samples (e.g.,
#' Sanger-sequenced clones or individual genomes). In deep
#' sequencing, reads are not independent lineages; they are
#' technical replicates of a viral population whose effective
#' size is much smaller than the read depth. Using raw read
#' depth (e.g., n = 2000) inflates the harmonic number
#' \eqn{a_1 = \sum 1/i}, making \eqn{\theta_W} artificially
#' small and biasing D toward positive values.
#'
#' We recommend capping \code{n} at a biologically meaningful
#' effective sample size. Lauring (2020 \emph{Annu Rev Virol})
#' suggests n = 100 as a conservative proxy. The high-level
#' \code{\link{calcNeutralityTests}} method applies this cap
#' automatically (\code{n = min(depth, 100)}).
#'
#' \strong{Interpretation warning:} Deep-sequencing D values are
#' \emph{not} directly comparable to classical population-genetic
#' D computed from independently sampled haplotypes. We recommend
#' reporting D as a \emph{descriptive statistic} of the site
#' frequency spectrum (e.g., "the majority of samples show
#' negative D, consistent with an excess of rare variants"),
#' rather than as a formal neutrality test with p-value
#' interpretation. For rigorous neutrality testing on deep
#' sequencing data, consider simulation-based approaches that
#' account for the read-sampling process.
#'
#' @references
#' Tajima F (1989). Statistical method for testing the neutral
#' mutation hypothesis by DNA polymorphism. \emph{Genetics}
#' 123:585-595.
#'
#' Lauring AS (2020). Within-Host Viral Diversity: A Window into
#' Viral Evolution. \emph{Annu Rev Virol} 7:63-81.
#' \doi{10.1146/annurev-virology-010320-061642}
#'
#' @export
#' @examples
#' tajimaD(S = 10, n = 1000, pi_hat = 0.001, genomeLength = 10000)
tajimaD <- function(S, n, pi_hat, genomeLength) {
    if (!is.numeric(S) || length(S) != 1L || S < 0)
        stop("'S' must be a single non-negative number.")
    if (!is.numeric(n) || length(n) != 1L || n <= 1)
        stop("'n' must be a single number > 1.")
    if (!is.numeric(pi_hat) || length(pi_hat) != 1L)
        stop("'pi_hat' must be a single numeric value.")
    if (!is.numeric(genomeLength) || length(genomeLength) != 1L ||
        genomeLength <= 0)
        stop("'genomeLength' must be a single positive number.")

    if (S < 3) {
        warning("Tajima's D needs at least three segregating sites; ",
                "with S = ", S, " it is returned as NA.", call. = FALSE)
        return(NA_real_)
    }

    n <- as.integer(round(n))
    if (n < 4L) {
        warning("Tajima's D needs a sample size of at least four; with ",
                "n = ", n, " its variance is zero and it is returned as NA.",
                call. = FALSE)
        return(NA_real_)
    }
    a1 <- sum(1 / seq_len(n - 1L))
    a2 <- sum(1 / (seq_len(n - 1L))^2)

    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))

    c1 <- b1 - 1 / a1
    c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2

    e1 <- c1 / a1
    e2 <- c2 / (a1^2 + a2)

    theta_W <- S / a1
    d <- pi_hat * genomeLength - theta_W
    denom <- sqrt(e1 * S + e2 * S * (S - 1))

    d / denom
}

#' Fu's Fs statistic
#'
#' Computes Fu's Fs (Fu 1997), a test of neutrality based on the
#' probability of observing at least the observed number of
#' alleles given the estimate of theta from nucleotide diversity.
#' Under the infinite-sites model:
#' \deqn{F_s = \ln\!\left(\frac{S'}{1 - S'}\right)}
#' where \eqn{S'} is the probability of observing a number of
#' alleles greater than or equal to the observed value, given
#' \eqn{\theta = \hat{\pi}}. A large negative Fs indicates an
#' excess of alleles (rare variants), consistent with population
#' expansion or genetic hitchhiking.
#'
#' The probability is approximated using Ewens' sampling formula
#' for the expected number of alleles.
#'
#' @param S Integer scalar. Number of segregating sites (alleles
#'   minus 1).
#' @param n Numeric scalar. Sample size. As in \code{\link{tajimaD}},
#'   do \strong{not} use raw read depth: use an effective sample size
#'   that reflects the number of independently sampled genomes.
#' @param theta_pi Numeric scalar. Theta estimated from
#'   nucleotide diversity (\eqn{\hat{\pi} \times L}).
#'
#' @return Numeric scalar. Fu's Fs, or \code{NA} with a warning when
#'   \eqn{S + 1 > n}: the statistic uses \eqn{S + 1} as the number of
#'   haplotypes, which cannot exceed the number of sequences sampled.
#'
#' @references
#' Fu YX (1997). Statistical tests of neutrality of mutations
#' against population growth, hitchhiking and background
#' selection. \emph{Genetics} 147:915-925.
#'
#' @export
#' @examples
#' fusFs(S = 10, n = 100, theta_pi = 5)
fusFs <- function(S, n, theta_pi) {
    if (!is.numeric(S) || length(S) != 1L || S < 0)
        stop("'S' must be a single non-negative number.")
    if (!is.numeric(n) || length(n) != 1L || n <= 1)
        stop("'n' must be a single number > 1.")
    if (!is.numeric(theta_pi) || length(theta_pi) != 1L ||
        theta_pi < 0)
        stop("'theta_pi' must be a single non-negative number.")

    if (S == 0 || theta_pi == 0) return(0)

    n <- as.integer(round(n))
    k_obs <- S + 1L
    if (k_obs > n) {
        warning("Fu's Fs is undefined when the number of segregating ",
                "sites (S = ", S, ") is not below the sample size (n = ",
                n, "): S + 1 haplotypes cannot occur among n sequences. ",
                "Returning NA.", call. = FALSE)
        return(NA_real_)
    }

    ## Compute P(K >= k_obs | theta) using Ewens sampling formula.
    ## Use recursive computation of Stirling numbers of the first
    ## kind (unsigned) for the exact probability.
    ## P(K = k | n, theta) = |s(n,k)| * theta^k / theta_(n)
    ## where theta_(n) = theta * (theta+1) * ... * (theta+n-1)

    ## Log of the rising factorial theta_(n)
    log_rising <- sum(log(theta_pi + seq_len(n) - 1))

    ## Compute log |s(n, k)| via recurrence:
    ## |s(n, k)| = (n-1) * |s(n-1, k)| + |s(n-1, k-1)|
    ## We only need k = k_obs .. n, so work on the full row.
    ## Use log-space to avoid overflow.

    ## For large n (deep sequencing), use the approximate Ewens
    ## expectation: E[K] = sum_{i=0}^{n-1} theta/(theta+i)
    ## and approximate the tail probability via a Gaussian.
    if (n > 500L) {
        ## Expected number of alleles under Ewens
        e_k <- sum(theta_pi / (theta_pi + seq_len(n) - 1))
        ## Variance: sum theta*i / (theta+i)^2
        v_k <- sum(theta_pi * (seq_len(n) - 1) /
                        (theta_pi + seq_len(n) - 1)^2)
        if (v_k <= 0) return(0)
        ## Continuity-corrected normal approximation
        z <- (k_obs - 0.5 - e_k) / sqrt(v_k)
        s_prime <- stats::pnorm(z, lower.tail = FALSE)
    } else {
        ## Exact computation via log Stirling numbers
        ## log_s[k] = log |s(n, k)| for k = 1..n
        ## Base case: |s(1,1)| = 1
        log_s <- rep(-Inf, n)
        log_s[1] <- 0  # log |s(1,1)| = 0
        for (i in 2:n) {
            new_log_s <- rep(-Inf, n)
            for (k in seq_len(i)) {
                ## |s(i,k)| = (i-1)*|s(i-1,k)| + |s(i-1,k-1)|
                term1 <- if (log_s[k] > -Inf) {
                    log(i - 1) + log_s[k]
                } else {
                    -Inf
                }
                term2 <- if (k >= 2 && log_s[k - 1] > -Inf) {
                    log_s[k - 1]
                } else {
                    -Inf
                }
                if (term1 == -Inf && term2 == -Inf) {
                    new_log_s[k] <- -Inf
                } else {
                    mx <- max(term1, term2)
                    new_log_s[k] <- mx + log(
                        exp(term1 - mx) + exp(term2 - mx))
                }
            }
            log_s <- new_log_s
        }
        ## P(K = k) = |s(n,k)| * theta^k / theta_(n)
        ## Sum P(K >= k_obs)
        ks <- seq(k_obs, n)
        log_probs <- vapply(ks, function(k) {
            log_s[k] + k * log(theta_pi) - log_rising
        }, numeric(1))
        max_lp <- max(log_probs)
        s_prime <- exp(max_lp) * sum(exp(log_probs - max_lp))
    }

    ## Clamp to avoid log(0) or log(Inf)
    s_prime <- max(s_prime, .Machine$double.xmin)
    s_prime <- min(s_prime, 1 - .Machine$double.eps)

    log(s_prime / (1 - s_prime))
}

# ============================================================
# High-level interface (operates on WHE)
# ============================================================

#' Compute neutrality test statistics per sample
#'
#' Calculates neutrality test statistics per sample from a
#' \code{\link{WithinHostExperiment}}. Implements Tajima's D
#' (Tajima 1989) and Fu's Fs (Fu 1997), adapted for deep
#' sequencing iSNV data.
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param indices Character vector. Which tests to compute:
#'   \code{"tajimaD"}, \code{"fusFs"} (default: both).
#' @param usePassedOnly Logical. If TRUE (default), only use
#'   variants where \code{qcPass == TRUE}.
#' @param ... Not used; supplying anything is an error.
#'
#' @details
#' No genome length is needed. Both statistics depend on the nucleotide
#' diversity summed over the genome, \eqn{\hat{\pi} L}, which is computed
#' directly from the frequencies, so a per-site value would only be
#' multiplied back by the length it was divided by.
#'
#' @return A \code{DataFrame} with one row per sample and columns
#'   \code{sample_id}, \code{n_sites} (the number of segregating sites
#'   each statistic rests on, that is, the variant sites present in the
#'   sample after QC) and one column per requested test. Tajima's D is
#'   \code{NA} for a sample with fewer than three segregating sites (see
#'   \code{\link{tajimaD}}); one warning names how many samples that
#'   affected.
#'
#' @references
#' Tajima F (1989). Statistical method for testing the neutral
#' mutation hypothesis by DNA polymorphism. \emph{Genetics}
#' 123:585-595.
#'
#' Fu YX (1997). Statistical tests of neutrality of mutations
#' against population growth, hitchhiking and background
#' selection. \emph{Genetics} 147:915-925.
#'
#' @export
#' @rdname calcNeutralityTests
#' @aliases calcNeutralityTests,WithinHostExperiment-method
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' calcNeutralityTests(whe)
setMethod("calcNeutralityTests", "WithinHostExperiment",
    function(x,
             indices = c("tajimaD", "fusFs"),
             usePassedOnly = TRUE, ...) {
    whe <- x
    .stop_on_unused_dots("calcNeutralityTests", ...)
    indices <- match.arg(indices, several.ok = TRUE)

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if ("qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    depth_mat <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids))
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    n_samples <- ncol(freq_mat)

    results <- lapply(seq_len(n_samples), function(j) {
        freq_j <- freq_mat[, j]
        if (usePassedOnly && !is.null(qc_mat))
            freq_j[!qc_mat[, j]] <- NA
        valid <- freq_j[!is.na(freq_j)]

        S <- length(valid)
        md <- if (!is.null(depth_mat)) {
            median(depth_mat[, j], na.rm = TRUE)
        } else {
            100
        }
        ## Cap at 100: deep sequencing reads are not independent
        ## lineages; using raw depth (e.g., 2000x) inflates a1
        ## and biases D toward positive values.
        n <- min(max(as.integer(round(md)), 2L), 100L)

        ## Nucleotide diversity summed over the genome (pi times L);
        ## both statistics use only this product, so no length is needed
        theta_pi <- if (S > 0) {
            2 * sum(valid * (1 - valid)) * n / (n - 1L)
        } else {
            0
        }

        vals <- list(sample_id = sample_ids[j], n_sites = S)

        if ("tajimaD" %in% indices) {
            ## below three sites tajimaD() warns per sample; one summary
            ## warning is given for the whole object instead
            vals$tajimaD <- if (S < 3L) NA_real_ else
                tajimaD(S, n, pi_hat = theta_pi, genomeLength = 1)
        }
        if ("fusFs" %in% indices)
            vals$fusFs <- fusFs(S, n, theta_pi)

        as.data.frame(vals, stringsAsFactors = FALSE)
    })

    out <- DataFrame(do.call(rbind, results))
    if ("tajimaD" %in% indices) {
        few <- out$n_sites < 3L
        if (any(few))
            warning("Tajima's D is NA for ", sum(few), " of ", nrow(out),
                    " sample(s), which have fewer than three segregating ",
                    "sites.", call. = FALSE)
    }
    out
})
