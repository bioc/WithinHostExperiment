# R/sfs.R
# Within-host Site Frequency Spectrum (SFS) as a first-class object
#
# Design rationale:
# The SFS is the fundamental summary statistic for population genetics.
# In within-host pathogen analysis, the SFS from deep sequencing has
# unique properties: (1) the "sample size" n is the read depth, not
# the number of individuals; (2) low-frequency variants are under-
# ascertained due to the detection threshold; (3) near-fixed variants
# (freq > 0.5) represent consensus-level changes, not iSNVs.
# This class makes all of these explicit and correctable.

#' @include AllClasses.R
#' @include methods.R
#' @importFrom methods setClass setValidity new is slot
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom stats median pnorm dbinom
NULL

# ============================================================
# WithinHostSFS class
# ============================================================

#' WithinHostSFS: Site Frequency Spectrum for within-host variants
#'
#' An S4 class representing the site frequency spectrum (SFS) of
#' within-host variants from deep sequencing data. The SFS bins
#' variant allele frequencies into discrete categories and stores
#' metadata needed for population-genetic inference.
#'
#' @slot counts Integer vector. Number of iSNVs in each frequency
#'   bin. Length equals \code{length(breaks) - 1}.
#' @slot breaks Numeric vector. Bin boundaries. Bins are
#'   right-closed, and the lowest break is included in the first
#'   bin. Default: 20 equal bins from the detection threshold to
#'   the maximum frequency.
#' @slot folded Logical scalar. If \code{TRUE}, the SFS is folded
#'   (minor allele frequency, 0 to 0.5); if \code{FALSE}, unfolded:
#'   the frequency of the alternative allele, 0 to 1.
#' @slot polarization Character scalar: which allele is treated as
#'   ancestral. \code{"none"} for a folded spectrum, which carries no
#'   polarisation; \code{"reference"} for an unfolded spectrum built by
#'   \code{\link{buildSFS}}, whose frequencies are those of the
#'   alternative allele, so the reference allele stands in for the
#'   ancestral one; \code{"outgroup"} for a spectrum whose frequencies
#'   are of alleles derived relative to an outgroup, which this package
#'   does not build but accepts. Statistics that need derived-allele
#'   frequencies, such as Fay and Wu's H, refuse \code{"none"}.
#' @slot sampleId Character scalar. Sample identifier.
#' @slot genomeLength Integer scalar. Reference genome length.
#' @slot nSites Integer scalar. Number of variant sites counted in
#'   the spectrum, that is, with a frequency between the lowest and the
#'   highest break; \code{\link{buildSFS}} counts only segregating sites,
#'   whose frequency lies strictly between 0 and 1.
#' @slot meanDepth Numeric scalar. Typical sequencing depth;
#'   \code{\link{buildSFS}} stores the median depth across sites.
#' @slot threshold Numeric scalar. Detection threshold applied
#'   (e.g., 0.03). Sites below this are unobservable.
#' @slot corrected Logical scalar. Whether ascertainment bias
#'   correction has been applied.
#' @slot maxFreq Numeric scalar, above 0 and at most 1 (default 1). The
#'   highest alternative allele frequency a call could have kept through
#'   QC: \code{\link{buildSFS}} stores the smallest \code{maxFreq} among
#'   the filters applied with \code{\link{flagISNV}} when it uses
#'   QC-passed calls only, and 1 otherwise. Below 1, the spectrum holds no
#'   call above it, so the high-frequency classes are empty by
#'   construction; \code{\link{neutralityFromSFS}} refuses such an
#'   unfolded spectrum, and \code{\link{compareSFS}} compares only spectra
#'   with the same value.
#'
#' @return An object of class \code{WithinHostSFS}, normally created by
#'   \code{\link{buildSFS}}; the \code{show} method prints a summary of the
#'   spectrum and returns \code{NULL} invisibly.
#'
#' @seealso \code{\link{buildSFS}} for constructing from a
#'   \code{WithinHostExperiment}, \code{\link{neutralityFromSFS}}
#'   for computing all neutrality tests from the SFS.
#'
#' @examples
#' data(example_whe)
#' ## a folded spectrum of the first sample, in ten bins from 3% to 50%;
#' ## the summary reports its polarisation, "none" for a folded spectrum,
#' ## and the genome length, read from seqlengths(example_whe)
#' sfs <- buildSFS(example_whe, sampleIdx = 1L, nBins = 10L)
#' sfs
#'
#' @aliases WithinHostSFS
#' @exportClass WithinHostSFS
.WithinHostSFS <- setClass("WithinHostSFS",
    slots = list(
        counts       = "integer",
        breaks       = "numeric",
        folded       = "logical",
        sampleId     = "character",
        genomeLength = "integer",
        nSites       = "integer",
        meanDepth    = "numeric",
        threshold    = "numeric",
        corrected    = "logical",
        polarization = "character",
        maxFreq      = "numeric"
    ),
    prototype = list(
        counts       = integer(0),
        breaks       = numeric(0),
        folded       = FALSE,
        sampleId     = NA_character_,
        genomeLength = NA_integer_,
        nSites       = 0L,
        meanDepth    = NA_real_,
        threshold    = 0.03,
        corrected    = FALSE,
        polarization = "reference",
        maxFreq      = 1
    )
)

setValidity("WithinHostSFS", function(object) {
    msg <- character()
    if (length(object@counts) != length(object@breaks) - 1L)
        msg <- c(msg, "length(counts) must equal length(breaks) - 1")
    if (any(object@counts < 0L))
        msg <- c(msg, "counts must be non-negative")
    if (is.unsorted(object@breaks))
        msg <- c(msg, "breaks must be sorted")
    pol <- object@polarization
    if (length(pol) != 1L || is.na(pol) ||
        !pol %in% c("none", "reference", "outgroup")) {
        msg <- c(msg, paste("polarization must be one of \"none\",",
                            "\"reference\" or \"outgroup\""))
    } else if (isTRUE(object@folded) && pol != "none") {
        msg <- c(msg, paste("a folded spectrum carries no polarisation:",
                            "polarization must be \"none\""))
    } else if (isFALSE(object@folded) && pol == "none") {
        msg <- c(msg, paste("an unfolded spectrum must say which allele is",
                            "ancestral: polarization \"reference\" or",
                            "\"outgroup\""))
    }
    mf <- object@maxFreq
    if (length(mf) != 1L || is.na(mf) || mf <= 0 || mf > 1)
        msg <- c(msg, "maxFreq must be a single number above 0 and at most 1")
    if (length(msg) == 0L) TRUE else msg
})

#' @rdname WithinHostSFS-class
#' @param object A \code{WithinHostSFS} object.
#' @export
setMethod("show", "WithinHostSFS", function(object) {
    cat("WithinHostSFS with", sum(object@counts), "variants in",
        length(object@counts), "bins\n")
    cat("  sample:", object@sampleId, "\n")
    cat("  folded:", object@folded, "\n")
    cat("  polarization:", object@polarization, "\n")
    cat("  threshold:", object@threshold, "\n")
    cat("  corrected:", object@corrected, "\n")
    if (!is.na(object@genomeLength))
        cat("  genome:", object@genomeLength, "bp\n")
    if (!is.na(object@meanDepth))
        cat("  median depth:", round(object@meanDepth), "x\n")
    if (object@maxFreq < 1)
        cat("  maxFreq:", object@maxFreq, "(QC removed every call above it)\n")
})


# ============================================================
# buildSFS -- construct from WithinHostExperiment
# ============================================================

## The one column a spectrum is built from: a sample named exactly once
## among the column names, or a position within them. Anything else,
## including several samples at once, is an error that lists the samples.
.sfsSampleColumn <- function(sampleIdx, samples, n) {
    if (length(sampleIdx) == 1L && !is.na(sampleIdx)) {
        if (is.character(sampleIdx)) {
            hit <- which(samples == sampleIdx)
            if (length(hit) == 1L) return(hit)
        } else if (is.numeric(sampleIdx) && sampleIdx == round(sampleIdx) &&
                   sampleIdx >= 1 && sampleIdx <= n) {
            return(as.integer(sampleIdx))
        }
    }
    shown <- if (length(sampleIdx) == 1L) paste0(" (got ", sampleIdx, ")")
             else paste0(" (got ", length(sampleIdx), " values)")
    stop("'sampleIdx' must give one sample, by a name that occurs once ",
         "among the column names or by a position from 1 to ", n, shown,
         ". Samples: ",
         if (length(samples)) paste(samples, collapse = ", ") else "unnamed",
         ".", call. = FALSE)
}

#' Build a site frequency spectrum from a WithinHostExperiment
#'
#' Constructs a \code{\link{WithinHostSFS}} object from the allele
#' frequencies in a \code{\link{WithinHostExperiment}}. One SFS is
#' built per sample. Only segregating sites enter it, those whose
#' frequency lies strictly between 0 and 1, as in
#' \code{\link{calcNeutralityTests}}: a frequency of 1, a fixed difference
#' from the reference, would otherwise fill the top bin of an unfolded
#' spectrum.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param sampleIdx Integer or character. The one sample (column) to
#'   use, by position or by a column name that occurs once. Default: 1.
#' @param nBins Integer. Number of frequency bins, a whole number of at
#'   least 1 (default: 20).
#' @param fold Logical. If \code{TRUE}, fold the SFS to minor
#'   allele frequency (0 to 0.5), giving \code{polarization = "none"}.
#'   If \code{FALSE}, bin the alternative allele frequency from 0 to 1,
#'   treating the reference allele as ancestral
#'   (\code{polarization = "reference"}). Default: \code{TRUE}.
#' @param usePassedOnly Logical. If \code{TRUE}, only use QC-passed
#'   variants. Default: \code{TRUE}. The spectrum's \code{maxFreq} slot
#'   then records the smallest \code{maxFreq} among the filters applied
#'   with \code{\link{flagISNV}}, above which no call can appear; with
#'   the default \code{\link{ISNVFilter}} that is 0.5, so an unfolded
#'   spectrum is empty above 0.5. Flag with \code{maxFreq = 1} to keep
#'   the high-frequency classes.
#' @param threshold Numeric. Detection threshold (default: 0.03).
#'   The first bin starts at this value, not 0. It must be at least 0
#'   and below the top of the spectrum: 0.5 when folded, 1 otherwise.
#' @param genomeLength Integer (optional). The genome length that the
#'   per-site statistics of \code{\link{neutralityFromSFS}} divide by. By
#'   default it is read from \code{seqlengths(rowRanges(whe))}, summed
#'   over every recorded sequence as in \code{\link{calcDiversity}}, and
#'   is \code{NA} when that is not recorded. The spectrum itself does not
#'   depend on it.
#'
#' @return A \code{\link{WithinHostSFS}} object.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L)
#' sfs
buildSFS <- function(whe, sampleIdx = 1L, nBins = 20L,
                     fold = TRUE, usePassedOnly = TRUE,
                     threshold = 0.03, genomeLength = NULL) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    freq_mat <- assay(whe, "altFreq")
    sampleIdx <- .sfsSampleColumn(sampleIdx, colnames(freq_mat),
                                  ncol(freq_mat))
    if (!is.numeric(nBins) || length(nBins) != 1L || !is.finite(nBins) ||
        nBins < 1 || nBins != round(nBins))
        stop("'nBins' must be a single whole number of at least 1.",
             call. = FALSE)
    top <- if (fold) 0.5 else 1
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0 || threshold >= top)
        stop("'threshold' must be a single number from 0 up to, but not ",
             "including, ", top, ", the top of ",
             if (fold) "a folded" else "an unfolded", " spectrum.",
             call. = FALSE)
    freq <- freq_mat[, sampleIdx]

    ## Apply QC, recording the frequency cap it imposed: the strictest
    ## maxFreq among the filters applied, or 1 when none was or QC is not
    ## being applied
    cap <- 1
    if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        qc <- assay(whe, "qcPass")[, sampleIdx]
        freq[!qc] <- NA
        filters <- qcFilters(whe)
        if (length(filters))
            cap <- min(1, vapply(filters, function(f)
                as.numeric(slot(f, "maxFreq")), numeric(1)))
    }
    ## only segregating sites, 0 < f < 1, as in calcNeutralityTests(): a
    ## fixed site would otherwise fill the top bin of an unfolded spectrum
    freq <- freq[!is.na(freq) & freq > 0 & freq < 1]

    ## Get sample metadata
    cd <- colData(whe)
    sid <- as.character(cd$sample_id[sampleIdx])

    ## Median depth across sites, stored in the meanDepth slot
    md <- if ("totalDepth" %in% assayNames(whe)) {
        median(assay(whe, "totalDepth")[, sampleIdx], na.rm = TRUE)
    } else {
        NA_real_
    }

    ## Fold if requested
    if (fold) {
        ## pmin() keeps an empty vector numeric, which ifelse() does not
        freq <- pmin(freq, 1 - freq)
        max_f <- 0.5
    } else {
        max_f <- 1.0
    }

    ## Build bins from threshold to max
    breaks <- seq(threshold, max_f, length.out = nBins + 1L)
    counts <- as.integer(hist(freq[freq >= threshold & freq <= max_f],
                              breaks = breaks, plot = FALSE)$counts)

    gl <- as.integer(if (!is.null(genomeLength)) genomeLength else
        .recordedGenomeLength(whe))

    new("WithinHostSFS",
        counts       = counts,
        breaks       = breaks,
        folded       = fold,
        sampleId     = sid,
        genomeLength = gl,
        nSites       = as.integer(sum(counts)),
        meanDepth    = md,
        threshold    = threshold,
        corrected    = FALSE,
        polarization = if (fold) "none" else "reference",
        maxFreq      = cap)
}


#' Build SFS for all samples in a WithinHostExperiment
#'
#' Convenience wrapper that calls \code{\link{buildSFS}} for each
#' sample and returns a named list of \code{\link{WithinHostSFS}}
#' objects.
#'
#' @inheritParams buildSFS
#'
#' @return A named list of \code{\link{WithinHostSFS}} objects.
#'
#' @export
#' @examples
#' data(example_whe)
#' sfs_list <- buildSFSList(example_whe, genomeLength = 1000L)
#' length(sfs_list)
buildSFSList <- function(whe, nBins = 20L, fold = TRUE,
                         usePassedOnly = TRUE, threshold = 0.03,
                         genomeLength = NULL) {
    n <- ncol(whe)
    sfs_list <- lapply(seq_len(n), function(j) {
        buildSFS(whe, sampleIdx = j, nBins = nBins, fold = fold,
                 usePassedOnly = usePassedOnly, threshold = threshold,
                 genomeLength = genomeLength)
    })
    names(sfs_list) <- as.character(colData(whe)$sample_id)
    sfs_list
}


# ============================================================
# Ascertainment bias correction
# ============================================================

#' Correct ascertainment bias in the within-host SFS
#'
#' Deep sequencing imposes a detection threshold (e.g., 3\%): sites
#' with true frequency below this threshold are unobservable. This
#' inflates the apparent fraction of moderate-frequency variants
#' relative to rare variants.
#'
#' This function applies a correction based on the expected number
#' of unobserved sites under a neutral beta distribution model:
#' given the observed SFS and depth-dependent detectability, it
#' estimates the missing mass in the lowest frequency bins.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object.
#' @param method Character. Correction method:
#'   \code{"truncation"} (default) rescales the whole spectrum for
#'   the unobservable part of the frequency range;
#'   \code{"binomial"} models the probability of detecting a
#'   variant at each frequency given the read depth.
#'
#' @return A corrected \code{\link{WithinHostSFS}} with
#'   \code{corrected = TRUE}. Counts are no longer raw variant
#'   counts, and \code{nSites} changes accordingly, which also
#'   changes the number of segregating sites used by
#'   \code{\link{neutralityFromSFS}}. Use uncorrected spectra for
#'   \code{\link{compareSFS}}.
#'
#' @details
#' The \code{"truncation"} method assumes a uniform density of
#' variant frequencies and multiplies every bin by
#' \eqn{f_{max} / (f_{max} - \mathrm{threshold})}, the reciprocal of
#' the observable fraction of the frequency range. It therefore
#' rescales the spectrum without changing its shape. The
#' \code{"binomial"} method uses the stored depth to compute
#' \eqn{P(\mathrm{detect} | \nu, n) =
#'   P(\mathrm{Alt} \ge k_{\min} | n, \nu)}{P(detect|freq,depth)}
#' where \eqn{k_{\min} =
#'   \lceil n \times \mathrm{threshold} \rceil}{k_min = ceil(n*threshold)},
#' and reweights each bin by the inverse of its detectability.
#' Detectability is floored at 0.01, so a bin is inflated at most
#' one hundredfold.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L)
#' sfs_corr <- correctSFSBias(sfs, method = "binomial")
#' sfs_corr
correctSFSBias <- function(sfs, method = c("truncation", "binomial")) {
    if (!is(sfs, "WithinHostSFS"))
        stop("'sfs' must be a WithinHostSFS object.")
    method <- match.arg(method)

    counts <- sfs@counts
    breaks <- sfs@breaks
    n_bins <- length(counts)
    thr    <- sfs@threshold
    depth  <- sfs@meanDepth

    if (method == "truncation") {
        ## Estimate the missing fraction: assume uniform density
        ## in the first bin extends below the threshold.
        ## fraction_observable = (max_freq - threshold) / max_freq
        max_f <- max(breaks)
        frac_obs <- (max_f - thr) / max_f
        ## Scale all counts by 1/frac_obs
        counts <- as.integer(round(counts / frac_obs))

    } else if (method == "binomial") {
        if (is.na(depth) || depth <= 0)
            stop("'binomial' correction requires mean depth ",
                 "in the SFS object.")
        n <- as.integer(round(depth))
        k_min <- ceiling(n * thr)

        ## For each bin midpoint, compute P(detect)
        mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
        p_detect <- vapply(mids, function(nu) {
            ## P(alt_reads >= k_min | n, nu)
            1 - pbinom(k_min - 1L, n, nu)
        }, numeric(1))

        ## Reweight: true_count = observed / P(detect)
        p_detect[p_detect < 0.01] <- 0.01  # floor to avoid Inf
        counts <- as.integer(round(counts / p_detect))
    }

    new("WithinHostSFS",
        counts       = counts,
        breaks       = breaks,
        folded       = sfs@folded,
        sampleId     = sfs@sampleId,
        genomeLength = sfs@genomeLength,
        nSites       = as.integer(sum(counts)),
        meanDepth    = sfs@meanDepth,
        threshold    = sfs@threshold,
        corrected    = TRUE,
        polarization = sfs@polarization,
        maxFreq      = sfs@maxFreq)
}


# ============================================================
# neutralityFromSFS -- all tests from one SFS
# ============================================================

#' Compute all neutrality statistics from a WithinHostSFS
#'
#' Computes Tajima's D, Fu's Fs, and Fay and Wu's H from a single
#' \code{\link{WithinHostSFS}} object. This is the recommended
#' interface for neutrality testing: build the SFS once, then
#' derive all statistics from it.
#'
#' \strong{Fay and Wu's H} measures an excess of high-frequency
#' derived variants, which can indicate positive selection or
#' genetic hitchhiking (Fay and Wu 2000). It needs the frequency of
#' the derived allele at each site, so the function stops on a folded
#' spectrum (\code{polarization = "none"}): folding maps a derived
#' allele at 0.9 and one at 0.1 to the same bin. For Tajima's D and
#' Fu's Fs alone, which need no polarisation, use
#' \code{\link{calcNeutralityTests}}. An unfolded spectrum built by
#' \code{\link{buildSFS}} has \code{polarization = "reference"}: its
#' frequencies are those of the alternative allele, so H here assumes
#' that the reference allele is ancestral at every site, which holds
#' only as far as the reference is close to the population the host was
#' infected from. H is defined as:
#' \deqn{H = \hat{\pi} - \hat{\theta}_H}
#' where \eqn{\hat{\pi} = \frac{2}{n(n-1)} \sum_i i (n - i) \xi_i}{pi =
#'   2/(n(n-1)) * sum(i * (n - i) * xi_i)} and \eqn{\hat{\theta}_H =
#'   \frac{2}{n(n-1)} \sum_i i^2 \xi_i}{theta_H = 2/(n(n-1)) * sum(i^2 * xi_i)}
#' are Fay and Wu's estimators on the same classes, and \eqn{\xi_i} is
#' the count of sites with \eqn{i} derived alleles. Each bin stands for
#' the class \eqn{i} nearest its midpoint, \code{round(mid * n)}, held
#' within 1 to \eqn{n - 1}. On matching classes both estimators have the
#' same expectation under neutrality, so H is centred on 0, and
#' \eqn{\hat{\pi}} carries the \eqn{n/(n-1)} correction that
#' \code{\link{calcNeutralityTests}} applies. A negative H indicates an
#' excess of high-frequency variants.
#'
#' Because H weights the high-frequency classes most, it is
#' uninterpretable when QC has removed them. The function therefore stops
#' on a spectrum whose \code{maxFreq} slot is below 1, as it is after
#' \code{\link{flagISNV}} with the default \code{\link{ISNVFilter}},
#' which drops every call above 0.5: H would then be positive whatever
#' the data. Flag the object with \code{ISNVFilter(maxFreq = 1)} before
#' \code{buildSFS(fold = FALSE)}.
#'
#' The per-site values need the spectrum's genome length, which
#' \code{\link{buildSFS}} takes from the object; a spectrum without one
#' is an error rather than a set of per-genome values labelled per site.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object.
#' @param n Integer. Effective sample size for neutrality tests.
#'   For deep sequencing data, do \strong{not} use raw depth.
#'   Default: \code{min(100, meanDepth)}, where the spectrum's
#'   \code{meanDepth} slot holds the median depth across sites, as
#'   recommended in the \code{\link{tajimaD}} documentation.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{S}{Number of segregating sites.}
#'     \item{n}{Effective sample size used.}
#'     \item{pi}{Nucleotide diversity per site, Fay and Wu's
#'       \eqn{\hat{\pi}} on the classes \eqn{i}; this is the value passed
#'       to \code{\link{tajimaD}} and, times the genome length, to
#'       \code{\link{fusFs}}.}
#'     \item{theta_W}{Watterson's theta per site.}
#'     \item{theta_H}{Fay and Wu's theta_H per site.}
#'     \item{tajimaD}{Tajima's D.}
#'     \item{fusFs}{Fu's Fs, with the number of haplotypes set to
#'       \eqn{S + 1}, the infinite-sites maximum. As \code{\link{fusFs}}
#'       explains, it then carries no haplotype information, depends only
#'       on \eqn{S}, \eqn{n} and \eqn{\hat{\pi}}, and lies well below 0
#'       under neutrality, so a negative value is not evidence of
#'       expansion.}
#'     \item{fayWuH}{Fay and Wu's H.}
#'     \item{sfs_corrected}{Whether the input SFS was bias-corrected.}
#'     \item{polarization}{The input's polarisation, \code{"reference"}
#'       or \code{"outgroup"}: which allele \code{theta_H} and
#'       \code{fayWuH} treat as ancestral.}
#'   }
#'
#' @references
#' Fay JC, Wu CI (2000). Hitchhiking under positive Darwinian
#' selection. \emph{Genetics} 155:1405-1413.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L, fold = FALSE)
#' neutralityFromSFS(sfs)
neutralityFromSFS <- function(sfs, n = NULL) {
    if (!is(sfs, "WithinHostSFS"))
        stop("'sfs' must be a WithinHostSFS object.")
    if (sfs@folded || identical(sfs@polarization, "none"))
        stop("Fay and Wu's H needs derived-allele frequencies, and this ",
             "spectrum is folded (polarization \"none\"). Build it with ",
             "buildSFS(fold = FALSE), or use calcNeutralityTests() for ",
             "Tajima's D and Fu's Fs, which need no polarisation.",
             call. = FALSE)
    if (sfs@maxFreq < 1)
        stop("Fay and Wu's H weights the high-frequency derived classes ",
             "most, and those were removed by QC: every call above ",
             "maxFreq = ", sfs@maxFreq, " was dropped before this spectrum ",
             "was built, so the classes are empty by construction and H ",
             "would be positive whatever the data. Flag the object with ",
             "flagISNV(x, ISNVFilter(..., maxFreq = 1)) and rebuild the ",
             "spectrum with buildSFS(fold = FALSE), or use ",
             "calcNeutralityTests() for Tajima's D and Fu's Fs.",
             call. = FALSE)

    S <- sfs@nSites
    if (is.null(n)) {
        n <- if (!is.na(sfs@meanDepth)) {
            min(100L, max(as.integer(round(sfs@meanDepth)), 2L))
        } else {
            100L
        }
    }
    n <- as.integer(n)
    if (n <= 1L) stop("'n' must be > 1.")

    gl <- sfs@genomeLength
    if (is.na(gl))
        stop("neutralityFromSFS() reports pi, Watterson's theta and ",
             "theta_H per site, which needs the genome length, and this ",
             "spectrum has none. Build it from an object whose ",
             "seqlengths(rowRanges()) are recorded, or pass genomeLength ",
             "to buildSFS().", call. = FALSE)

    ## Each bin stands for the derived-allele class i = round(mid * n),
    ## held within the segregating classes 1..n-1, and pi and theta_H are
    ## Fay and Wu's estimators on those same classes:
    ## theta_pi = sum 2 i (n - i) xi_i / (n (n - 1)), which carries the
    ## n/(n-1) correction, and theta_H = sum 2 i^2 xi_i / (n (n - 1)), so
    ## that H = theta_pi - theta_H has expectation 0 under neutrality
    breaks <- sfs@breaks
    mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
    counts <- sfs@counts
    cls <- pmin(pmax(round(mids * n), 1), n - 1)
    norm <- 2 / (as.numeric(n) * (n - 1))

    ## Nucleotide diversity from SFS
    pi_per_site <- if (S > 0 && gl > 0) {
        norm * sum(counts * cls * (n - cls)) / gl
    } else {
        0
    }

    ## Watterson's theta
    a1 <- sum(1 / seq_len(n - 1L))
    theta_W <- if (S > 0) (S / a1) / gl else 0

    ## Fay and Wu's theta_H
    theta_H <- if (S > 0) {
        norm * sum(counts * cls^2) / gl
    } else {
        0
    }

    ## Tajima's D
    D <- tajimaD(S, n, pi_per_site, gl)

    ## Fu's Fs
    theta_pi <- pi_per_site * gl
    Fs <- fusFs(S, n, theta_pi)

    ## Fay and Wu's H
    H <- pi_per_site - theta_H

    list(
        S             = S,
        n             = n,
        pi            = pi_per_site,
        theta_W       = theta_W,
        theta_H       = theta_H,
        tajimaD       = D,
        fusFs         = Fs,
        fayWuH        = H,
        sfs_corrected = sfs@corrected,
        polarization  = sfs@polarization
    )
}


# ============================================================
# compareSFS -- compare two SFS objects
# ============================================================

#' Compare two within-host SFS objects
#'
#' Computes summary statistics and a chi-squared test for
#' homogeneity between two \code{\link{WithinHostSFS}} objects
#' (e.g., before vs. after QC, or two timepoints).
#'
#' The test compares raw variant counts in matching bins, so the two
#' spectra must share their breaks, their folding, their polarisation
#' and their QC frequency cap (the \code{maxFreq} slot), and neither may
#' be bias-corrected: the counts of a
#' corrected spectrum are rescaled and no longer multinomial. The
#' function stops, naming the difference, when any of these does not
#' hold. Bins where both spectra are
#' empty are dropped, and when only two bins remain
#' \code{\link[stats]{chisq.test}} applies Yates' continuity
#' correction. No test is run, and the statistic, p-value and
#' diagnostics are all \code{NA}, when either spectrum is empty or fewer
#' than two bins hold variants. The chi-squared approximation is unreliable
#' when an expected count falls below 5; instead of a warning, that
#' condition is reported in \code{chisq_approx_ok}, and the p-value
#' should then be read with caution.
#'
#' @param sfs1 A \code{\link{WithinHostSFS}} object.
#' @param sfs2 A \code{\link{WithinHostSFS}} object with the same
#'   breaks, folding, polarisation and \code{maxFreq} as \code{sfs1}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{chisq_stat}{Chi-squared statistic, \code{NA} when no test
#'       is run.}
#'     \item{chisq_p}{Chi-squared p-value, \code{NA} when no test is
#'       run.}
#'     \item{df}{Degrees of freedom, \code{NA} when no test is run.}
#'     \item{min_expected}{Smallest expected count under homogeneity,
#'       or \code{NA} when either spectrum is empty or fewer than two
#'       bins remain, and no test is run.}
#'     \item{chisq_approx_ok}{\code{TRUE} when every expected count
#'       is at least 5, the condition under which
#'       \code{\link[stats]{chisq.test}} does not warn; \code{NA}
#'       alongside the other \code{NA} diagnostics when no test is
#'       run.}
#'     \item{n1}{Total variants in sfs1.}
#'     \item{n2}{Total variants in sfs2.}
#'     \item{proportions1}{Proportion in each bin, sfs1.}
#'     \item{proportions2}{Proportion in each bin, sfs2.}
#'   }
#'
#' @export
#' @examples
#' data(example_whe)
#' whe1 <- example_whe
#' sfs1 <- buildSFS(whe1, genomeLength = 1000L)
#' sfs2 <- buildSFS(whe1, sampleIdx = 2L, genomeLength = 1000L)
#' compareSFS(sfs1, sfs2)
compareSFS <- function(sfs1, sfs2) {
    if (!is(sfs1, "WithinHostSFS") || !is(sfs2, "WithinHostSFS"))
        stop("Both arguments must be WithinHostSFS objects.")
    if (length(sfs1@counts) != length(sfs2@counts))
        stop("SFS objects must have the same number of bins.")
    if (!isTRUE(all.equal(sfs1@breaks, sfs2@breaks)))
        stop("SFS objects must share their bin breaks; the counts of ",
             "differently placed bins cannot be compared.\n  sfs1: ",
             paste(signif(range(sfs1@breaks), 3), collapse = " to "),
             "\n  sfs2: ",
             paste(signif(range(sfs2@breaks), 3), collapse = " to "),
             call. = FALSE)
    if (!identical(sfs1@folded, sfs2@folded) ||
        !identical(sfs1@polarization, sfs2@polarization))
        stop("SFS objects must both be folded or both unfolded with the ",
             "same polarisation (sfs1: ",
             if (sfs1@folded) "folded" else sfs1@polarization, "; sfs2: ",
             if (sfs2@folded) "folded" else sfs2@polarization, ").",
             call. = FALSE)
    if (!isTRUE(all.equal(sfs1@maxFreq, sfs2@maxFreq)))
        stop("SFS objects must share their QC frequency cap; calls above ",
             "maxFreq were removed from one spectrum but not from the ",
             "other (sfs1: ", sfs1@maxFreq, "; sfs2: ", sfs2@maxFreq, ").",
             call. = FALSE)
    if (sfs1@corrected || sfs2@corrected)
        stop("compareSFS() tests raw counts; a bias-corrected spectrum ",
             "holds rescaled counts. Compare the spectra before ",
             "correctSFSBias().", call. = FALSE)

    c1 <- sfs1@counts
    c2 <- sfs2@counts
    n1 <- sum(c1)
    n2 <- sum(c2)

    ## Remove bins where both are zero. An empty spectrum carries no
    ## shape to compare, and its row of expected counts would be all zero,
    ## so it takes the same no-test path as fewer than two bins.
    keep <- (c1 + c2) > 0
    if (n1 == 0L || n2 == 0L || sum(keep) < 2L) {
        return(list(chisq_stat = NA_real_, chisq_p = NA_real_,
                    df = NA_integer_, min_expected = NA_real_,
                    chisq_approx_ok = NA, n1 = n1, n2 = n2,
                    proportions1 = c1 / max(n1, 1),
                    proportions2 = c2 / max(n2, 1)))
    }

    mat <- rbind(c1[keep], c2[keep])
    min_expected <- min(outer(rowSums(mat), colSums(mat)) / sum(mat))
    ## The small-expected-count warning is returned as chisq_approx_ok;
    ## any other warning from chisq.test() still reaches the caller.
    approx_warning <- gettext(
        "Chi-squared approximation may be incorrect", domain = "R-stats"
    )
    test <- withCallingHandlers(
        chisq.test(mat),
        warning = function(w) {
            if (identical(conditionMessage(w), approx_warning))
                invokeRestart("muffleWarning")
        }
    )

    list(
        chisq_stat      = unname(test$statistic),
        chisq_p         = test$p.value,
        df              = unname(test$parameter),
        min_expected    = min_expected,
        chisq_approx_ok = min_expected >= 5,
        n1              = n1,
        n2              = n2,
        proportions1    = c1 / max(n1, 1),
        proportions2    = c2 / max(n2, 1)
    )
}


# ============================================================
# plotSFS -- publication-quality SFS visualisation
# ============================================================

#' Plot a within-host site frequency spectrum
#'
#' Creates a publication-quality bar plot of a
#' \code{\link{WithinHostSFS}} object.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object, or a named
#'   list of \code{WithinHostSFS} objects for overlay.
#' @param normalize Logical. If \code{TRUE}, plot proportions
#'   instead of counts. Default: \code{FALSE}.
#'
#' @return A ggplot object.
#'
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     vcf <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(vcf,
#'         colData = S4Vectors::DataFrame(sample_id = "donor"),
#'         caller = "ivar")
#'     sfs <- buildSFS(whe, genomeLength = 1000L)
#'     plotSFS(sfs)
#' }
plotSFS <- function(sfs, normalize = FALSE) {
    .requireGgplot2()

    if (is(sfs, "WithinHostSFS")) {
        sfs <- list(sfs)
        names(sfs) <- sfs[[1]]@sampleId
    }
    if (!is.list(sfs) || !all(vapply(sfs, is, logical(1), "WithinHostSFS")))
        stop("'sfs' must be a WithinHostSFS or a named list of them.")

    df_list <- lapply(names(sfs), function(nm) {
        s <- sfs[[nm]]
        breaks <- s@breaks
        mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
        y <- if (normalize) s@counts / max(sum(s@counts), 1) else s@counts
        data.frame(
            freq = mids * 100,
            count = y,
            sample = nm,
            stringsAsFactors = FALSE)
    })
    plot_df <- do.call(rbind, df_list)
    n_samples <- length(sfs)

    ylab <- if (normalize) "Proportion" else "Number of iSNVs"

    p <- ggplot2::ggplot(plot_df,
        ggplot2::aes(x = .data$freq, y = .data$count))

    if (n_samples == 1) {
        p <- p + ggplot2::geom_col(fill = .whe_pal$blue,
            colour = "white", linewidth = 0.2, alpha = 0.90,
            width = diff(sfs[[1]]@breaks)[1] * 100 * 0.85)
    } else {
        p <- p + ggplot2::geom_col(
            ggplot2::aes(fill = .data$sample),
            colour = "white", linewidth = 0.2, alpha = 0.70,
            position = "dodge",
            width = diff(sfs[[1]]@breaks)[1] * 100 * 0.85) +
            ggplot2::scale_fill_manual(name = NULL,
                values = rep(unlist(.whe_pal, use.names = FALSE),
                             length.out = n_samples))
    }

    p <- p +
        ggplot2::scale_x_continuous(
            ifelse(sfs[[1]]@folded,
                   "Minor allele frequency (%)",
                   "Alternative allele frequency (%)")) +
        ggplot2::labs(y = ylab) +
        .theme_pub(base_size = 9)

    p
}
