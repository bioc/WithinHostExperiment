# R/diversity.R
# Within-host diversity indices

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqlengths seqlevels
#' @importFrom stats median wilcox.test kruskal.test
NULL

## The genome length a per-site index divides by: the summed lengths of
## every sequence recorded for the reference in seqlengths(rowRanges(x)),
## whether or not a site lies on it, so that a segment without calls
## still counts as sequence in which none was found. A length that is not
## recorded is an error that says how to record it.
.genomeLengthOf <- function(x, indices) {
    rr <- rowRanges(x)
    lv <- seqlevels(rr)
    lens <- seqlengths(rr)[lv]
    missing_len <- lv[is.na(lens)]
    if (length(lv) == 0L || length(missing_len)) {
        named <- if (length(missing_len)) missing_len else "<seqname>"
        stop("The index ", paste(sQuote(indices, FALSE), collapse = ", "),
             " divides by the genome length, which is not recorded",
             if (length(missing_len)) paste0(" for ",
                 paste(missing_len, collapse = ", ")),
             " in seqlengths(rowRanges(x)). Record it from the reference ",
             "the calls were made against: read the files with ",
             "readWithinHost(genome = <DNAStringSet>), or set it with ",
             "seqlengths(x)[c(", paste0('"', named, '"', collapse = ", "),
             ")] <- <length>.",
             if (length(missing_len))
                 paste0(" A sequence that is not part of this genome can ",
                        "be dropped with x <- GenomeInfoDb::keepSeqlevels(",
                        "x, GenomeInfoDb::seqlevelsInUse(x), pruning.mode = ",
                        "\"coarse\") if it carries no call."),
             call. = FALSE)
    }
    sum(as.numeric(lens))
}

## The same length, or NA when it is not fully recorded; for objects such
## as a spectrum that can be built without it.
.recordedGenomeLength <- function(x) {
    lens <- seqlengths(rowRanges(x))
    if (length(lens) == 0L || anyNA(lens)) return(NA_real_)
    sum(as.numeric(lens))
}

## The lowest frequency a call could have kept through QC: the strictest
## minFreq of the filters applied, or NA when none was, or when QC is
## not being applied.
.detectionLimitOf <- function(x, usePassedOnly) {
    filters <- qcFilters(x)
    if (!usePassedOnly || length(filters) == 0L) return(NA_real_)
    max(vapply(filters, function(f) as.numeric(slot(f, "minFreq")),
               numeric(1)))
}

## Allele frequencies for the index functions: numbers in [0, 1], with NA
## allowed (an all-NA logical vector counts as numeric). Text is an error
## rather than being coerced to NA, and so is a value outside the range.
.checkFreq <- function(freq) {
    if (!is.numeric(freq) && !(is.logical(freq) && all(is.na(freq))))
        stop("'freq' must be a numeric vector of allele frequencies.",
             call. = FALSE)
    bad <- !is.na(freq) & (freq < 0 | freq > 1)
    if (any(bad))
        stop("'freq' must lie between 0 and 1; ", sum(bad),
             " value(s) do not, the first being ", freq[bad][1L], ".",
             call. = FALSE)
    as.numeric(freq)
}

# ============================================================
# Standalone math functions (exported, work on plain vectors)
# ============================================================

#' Shannon entropy from iSNV frequencies
#'
#' Calculates per-site biallelic Shannon entropy, summed across
#' all iSNV sites:
#' \deqn{H = -\sum_i [p_i \log(p_i) + (1-p_i) \log(1-p_i)]}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}. This formulation treats each iSNV as an independent
#' biallelic locus, following McCrone & Lauring (2016) and Popa
#' et al. (2020). Higher values indicate greater within-host
#' diversity.
#'
#' @param freq Numeric vector. Alternative allele frequencies
#'   between 0 and 1; anything else is an error. Values of 0 and 1,
#'   which carry no diversity, and NA are excluded.
#'
#' @return Numeric scalar. Shannon entropy.
#'
#' @references
#' Lauring AS (2020). Within-Host Viral Diversity: A Window into
#' Viral Evolution. \emph{Annu Rev Virol} 7:63-81.
#' \doi{10.1146/annurev-virology-010320-061642}
#'
#' @export
#' @examples
#' shannonISNV(c(0.05, 0.15, 0.30, 0.45))
#' shannonISNV(0.5) # log(2)
shannonISNV <- function(freq) {
    freq <- .checkFreq(freq)
    freq <- freq[!is.na(freq) & freq > 0 & freq < 1]
    if (length(freq) == 0L) return(0)
    -sum(freq * log(freq) + (1 - freq) * log(1 - freq))
}

#' Nucleotide diversity (pi) from iSNV frequencies
#'
#' Estimates nucleotide diversity following Nei & Li (1979),
#' adapted for deep sequencing allele frequencies:
#' \deqn{\pi = \frac{n}{n-1} \cdot \frac{2}{L} \sum_i p_i(1 - p_i)}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}, \eqn{L} is the genome length, and \eqn{n} is the
#' number of reads sampled at a site, for which mean read depth
#' stands in. The \eqn{n/(n-1)} correction removes the
#' finite-sample bias of estimating a frequency from \eqn{n}
#' draws (Nei 1987, Eq. 10.5). This is a different use of
#' "sample size" from the one in \code{\link{tajimaD}}, where
#' reads must not stand in for independent lineages. When
#' \code{meanDepth} is not provided, the uncorrected estimator
#' is returned (equivalent to infinite sample size).
#'
#' @param freq Numeric vector. Alternative allele frequencies, between 0
#'   and 1; anything else is an error. NA values are dropped.
#' @param genomeLength Integer scalar. Reference genome length in bp.
#' @param meanDepth Numeric scalar (optional). Mean read depth for
#'   finite-sample correction. If \code{NULL}, no correction is
#'   applied.
#'
#' @return Numeric scalar. Nucleotide diversity.
#'
#' @references
#' Nei M, Li WH (1979). Mathematical model for studying genetic
#' variation in terms of restriction endonucleases. \emph{Proc
#' Natl Acad Sci USA} 76:5269-5273.
#' \doi{10.1073/pnas.76.10.5269}
#'
#' Nei M (1987). \emph{Molecular Evolutionary Genetics.}
#' Columbia University Press, New York. Chapter 10.
#'
#' @export
#' @examples
#' piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588)
#' piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588, meanDepth = 1000)
piISNV <- function(freq, genomeLength, meanDepth = NULL) {
    freq <- .checkFreq(freq)
    if (!is.numeric(genomeLength) || length(genomeLength) != 1L ||
        genomeLength <= 0)
        stop("'genomeLength' must be a single positive number.")
    freq <- freq[!is.na(freq)]
    if (length(freq) == 0L) return(0)
    raw_pi <- (2 / genomeLength) * sum(freq * (1 - freq))
    ## Apply finite-sample correction n/(n-1) (Nei 1987)
    if (!is.null(meanDepth)) {
        if (!is.numeric(meanDepth) || length(meanDepth) != 1L ||
            meanDepth <= 1)
            stop("'meanDepth' must be a single number > 1.")
        n <- as.integer(round(meanDepth))
        if (n < 2L)
            stop("'meanDepth' of ", meanDepth, " rounds to ", n,
                 "; the n/(n-1) correction needs a depth that rounds to ",
                 "2 or more.")
        raw_pi <- raw_pi * n / (n - 1L)
    }
    raw_pi
}

#' Watterson's theta from iSNV count
#'
#' Estimates Watterson's theta (Watterson 1975):
#' \deqn{\hat{\theta}_W = S / (a_n \cdot L)}
#' where \eqn{S} is the number of segregating sites,
#' \eqn{a_n = \sum_{i=1}^{n-1} 1/i} is the harmonic number,
#' and \eqn{L} is the genome length, with \eqn{n} the number of
#' sampled viral genomes.
#'
#' As \code{\link{tajimaD}} sets out, raw read depth should
#' \strong{not} be passed as \eqn{n} for deep sequencing data:
#' reads are technical replicates of a much smaller viral
#' population, and a depth-sized \eqn{n} inflates \eqn{a_n} and
#' makes \eqn{\theta_W} correspondingly small. Pass an effective
#' sample size, or a depth capped at one, as
#' \code{\link{calcNeutralityTests}} does.
#'
#' @param nSites Integer scalar. Number of segregating (iSNV) sites.
#' @param genomeLength Integer scalar. Genome length in bp.
#' @param meanDepth Numeric scalar. Sample size \eqn{n}: an
#'   effective number of sampled genomes, not raw read depth. See
#'   Details.
#'
#' @return Numeric scalar. Watterson's theta estimate.
#'
#' @references
#' Watterson GA (1975). On the number of segregating sites in
#' genetical models without recombination. \emph{Theor Popul
#' Biol} 7:256-276.
#' \doi{10.1016/0040-5809(75)90020-9}
#'
#' @export
#' @examples
#' ## n is an effective sample size, not the raw depth
#' wattersonISNV(15, genomeLength = 13588, meanDepth = 100)
wattersonISNV <- function(nSites, genomeLength, meanDepth) {
    if (!is.numeric(nSites) || length(nSites) != 1L || nSites < 0)
        stop("'nSites' must be a single non-negative number.")
    if (!is.numeric(genomeLength) || length(genomeLength) != 1L ||
        genomeLength <= 0)
        stop("'genomeLength' must be a single positive number.")
    if (!is.numeric(meanDepth) || length(meanDepth) != 1L ||
        meanDepth <= 1)
        stop("'meanDepth' must be a single number > 1.")
    if (nSites == 0) return(0)
    n <- max(as.integer(round(meanDepth)), 2L)
    a_n <- sum(1 / seq_len(n - 1L))
    nSites / (a_n * genomeLength)
}

#' Simpson diversity index for iSNV frequencies
#'
#' Calculates Simpson diversity summed over biallelic iSNV sites:
#' \deqn{D = \sum_i \left[1 - p_i^2 - (1 - p_i)^2\right] =
#'   \sum_i 2 p_i (1 - p_i)}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}. Each site contributes the probability that two alleles
#' drawn at random at that site differ (at most 0.5), so the index is
#' the expected number of differing sites. Higher values indicate
#' greater diversity.
#'
#' @param freq Numeric vector. Alternative allele frequencies
#'   between 0 and 1; anything else is an error. Values of 0 and 1,
#'   which carry no diversity, and NA are excluded.
#'
#' @return Numeric scalar. Simpson diversity index.
#'
#' @export
#' @examples
#' simpsonISNV(c(0.05, 0.15, 0.30, 0.45))
#' simpsonISNV(0.5) # maximum per-site diversity = 0.5
simpsonISNV <- function(freq) {
    freq <- .checkFreq(freq)
    freq <- freq[!is.na(freq) & freq > 0 & freq < 1]
    if (length(freq) == 0L) return(0)
    sum(1 - freq^2 - (1 - freq)^2)
}

#' Chao1 richness estimator for iSNV data
#'
#' Estimates total species (variant) richness using the Chao1
#' estimator (Chao 1984):
#' \deqn{\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2 f_2}}
#' where \eqn{S_{obs}} is the observed number of variants,
#' \eqn{f_1} is the number of singletons (variants with
#' \code{altCount == 1}), and \eqn{f_2} is the number of
#' doubletons (variants with \code{altCount == 2}). When
#' \eqn{f_2 = 0}, the bias-corrected form is used:
#' \eqn{S_{obs} + f_1 (f_1 - 1) / 2}.
#'
#' @param counts Integer vector. Alternative allele read counts
#'   (\code{altCount}) for each variant.
#' @param detected Integer scalar. Number of observed variants
#'   (\eqn{S_{obs}}).
#'
#' @return Numeric scalar. Chao1 richness estimate.
#'
#' @references
#' Chao A (1984). Nonparametric estimation of the number of
#' classes in a population. \emph{Scand J Statist} 11:265-270.
#'
#' @export
#' @examples
#' chao1ISNV(c(1, 1, 2, 5, 10), detected = 5)
chao1ISNV <- function(counts, detected) {
    counts <- as.integer(counts)
    counts <- counts[!is.na(counts)]
    if (length(counts) == 0L) return(as.numeric(detected))
    f1 <- sum(counts == 1L)
    f2 <- sum(counts == 2L)
    if (f2 == 0L) {
        detected + f1 * (f1 - 1) / 2
    } else {
        detected + f1^2 / (2 * f2)
    }
}

# ============================================================
# High-level interface (operates on WHE)
# ============================================================

#' Compute within-host diversity indices
#'
#' Calculates diversity metrics per sample from a
#' \code{\link{WithinHostExperiment}}. Implements standard
#' population genetics indices adapted for deep sequencing
#' iSNV data, following Lauring (2020).
#'
#' @references
#' Lauring AS (2020). Within-Host Viral Diversity: A Window into
#' Viral Evolution. \emph{Annu Rev Virol} 7:63-81.
#' \doi{10.1146/annurev-virology-010320-061642}
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param indices Character vector. Which indices to compute:
#'   \code{"shannon"}, \code{"pi"}, \code{"watterson"},
#'   \code{"richness"}, \code{"simpson"}, \code{"chao1"}
#'   (default: all six).
#' @param usePassedOnly Logical. If TRUE (default), only use variants
#'   where \code{qcPass == TRUE}.
#' @param ... Not used; supplying anything is an error.
#'
#' @details
#' \strong{Genome length.} \code{"pi"} and \code{"watterson"} are per
#' site and divide by the genome length, which is read from
#' \code{seqlengths(rowRanges(x))} as the summed length of every
#' sequence recorded there, including segments of a segmented genome on
#' which no site lies. It is a property of the reference the
#' calls were made against, not a number to type: \code{\link{readWithinHost}}
#' records it from the VCF header's \code{##contig} lines or from a
#' \code{DNAStringSet} passed as \code{genome}, and
#' \code{\link{annotateCodonChange}} from its reference. When it is
#' not recorded the function stops and says how to record it, for
#' example \code{seqlengths(x) <- c(seg = 29903L)}.
#'
#' For \code{"watterson"}, the sample's median read depth is passed as
#' the sample size \eqn{n}, which for the reasons given in
#' \code{\link{tajimaD}} makes \eqn{\theta_W} small;
#' \code{\link{calcNeutralityTests}} caps \eqn{n} instead, so the two
#' need not agree.
#'
#' \strong{Reading a zero.} A diversity of zero means something
#' different when the sample had no call at all, when every call failed
#' QC, and when calls were present but below the frequency QC allows.
#' Two columns keep these apart: \code{n_sites_assessed}, the number of
#' sites at which the sample has a frequency, before QC, and
#' \code{detection_limit}, the lowest frequency a call could have kept
#' through the filters applied (the strictest \code{minFreq} among
#' them), \code{NA} when no filter was applied or
#' \code{usePassedOnly = FALSE}. \code{richness} counts the sites that
#' entered the indices after QC: the segregating sites, whose frequency
#' lies strictly between 0 and 1. A frequency of 0, or of 1, a fixed
#' difference from the reference such as iVar reports, is not
#' segregating within the host, so it enters none of the indices; the
#' site count of \code{"watterson"} and the singletons and doubletons of
#' \code{"chao1"} come from the same sites. A site counts as assessed
#' only when the
#' sample has a call there: a VCF lists only the sites a caller reported,
#' so a site sequenced without a call cannot be told apart from one not
#' sequenced at all, and neither is counted.
#'
#' @return A \code{DataFrame} with one row per sample and columns
#'   \code{sample_id}, \code{n_sites_assessed}, \code{detection_limit}
#'   and one column per requested index.
#'
#' @export
#' @rdname calcDiversity
#' @aliases calcDiversity,WithinHostExperiment-method
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' ## the genome length comes from the reference, recorded on import
#' ref <- Biostrings::readDNAStringSet(system.file("extdata",
#'     "test_ref.fa", package = "WithinHostExperiment"))
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar", genome = ref)
#' calcDiversity(whe)
setMethod("calcDiversity", "WithinHostExperiment",
    function(x,
             indices = c("shannon", "pi",
                         "watterson", "richness",
                         "simpson", "chao1"),
             usePassedOnly = TRUE, ...) {
    whe <- x
    .stop_on_unused_dots("calcDiversity", ...)
    indices <- match.arg(indices, several.ok = TRUE)
    if (!is.logical(usePassedOnly) || length(usePassedOnly) != 1L ||
        is.na(usePassedOnly))
        stop("'usePassedOnly' must be TRUE or FALSE.")

    per_site <- intersect(c("pi", "watterson"), indices)
    genomeLength <- if (length(per_site)) .genomeLengthOf(whe, per_site)
                    else NULL
    detection_limit <- .detectionLimitOf(whe, usePassedOnly)

    freq_mat  <- assay(whe, "altFreq")
    qc_mat    <- if ("qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    if ("watterson" %in% indices &&
        !"totalDepth" %in% assayNames(whe)) {
        warning("No 'totalDepth' assay available; using default ",
                "meanDepth = 1000 for Watterson's theta. ",
                "This may produce inaccurate estimates.",
                call. = FALSE)
    }
    depth_mat <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")
    } else {
        NULL
    }
    if ("chao1" %in% indices &&
        !"altCount" %in% assayNames(whe)) {
        warning("No 'altCount' assay available; 'chao1' will ",
                "equal observed richness.",
                call. = FALSE)
    }
    alt_count_mat <- if ("altCount" %in% assayNames(whe)) {
        assay(whe, "altCount")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }
    n_samples <- ncol(freq_mat)

    results <- lapply(seq_len(n_samples), function(j) {
        freq_j <- freq_mat[, j]
        n_assessed <- sum(!is.na(freq_j))
        ## Apply QC mask
        if (usePassedOnly && !is.null(qc_mat)) {
            freq_j[!qc_mat[, j]] <- NA
        }
        ## the segregating sites: both alleles present, 0 < f < 1, as in
        ## shannonISNV() and simpsonISNV()
        seg <- !is.na(freq_j) & freq_j > 0 & freq_j < 1
        valid <- freq_j[seg]

        vals <- list(sample_id = sample_ids[j],
                     n_sites_assessed = n_assessed,
                     detection_limit = detection_limit)
        ## Median depth for this sample (used by pi and watterson)
        md <- if (!is.null(depth_mat)) {
            median(depth_mat[, j], na.rm = TRUE)
        } else {
            NULL
        }
        ## a sample with no depth anywhere, such as one with no calls, is
        ## treated as if there were no depth assay
        if (!is.null(md) && is.na(md)) md <- NULL

        if ("richness" %in% indices) vals$richness <- length(valid)
        if ("shannon" %in% indices) vals$shannon <- shannonISNV(valid)
        if ("pi" %in% indices) {
            vals$pi <- piISNV(valid, genomeLength,
                              meanDepth = md)
        }
        if ("watterson" %in% indices) {
            wd <- if (!is.null(md)) md else 1000
            vals$watterson <- wattersonISNV(length(valid),
                                            genomeLength, wd)
        }
        if ("simpson" %in% indices) vals$simpson <- simpsonISNV(valid)
        if ("chao1" %in% indices) {
            n_obs <- length(valid)
            if (!is.null(alt_count_mat)) {
                ac_j <- alt_count_mat[, j]
                if (usePassedOnly && !is.null(qc_mat)) {
                    ac_j[!qc_mat[, j]] <- NA
                }
                ac_valid <- ac_j[!is.na(ac_j) & seg]
                vals$chao1 <- chao1ISNV(ac_valid, n_obs)
            } else {
                vals$chao1 <- as.numeric(n_obs)
            }
        }
        as.data.frame(vals, stringsAsFactors = FALSE)
    })

    DataFrame(do.call(rbind, results))
})

#' Compare diversity between sample groups
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param groupBy Character. Column name in \code{colData} to group by.
#' @param index Character. Which index to compare (default: "shannon").
#'   \code{"pi"} and \code{"watterson"} need the genome length recorded
#'   in \code{seqlengths(rowRanges(whe))}; see
#'   \code{\link{calcDiversity}}.
#'
#' @return A \code{DataFrame} with test results.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' ## Single-group, returns diversity table:
#' compareDiversity(whe, groupBy = "sample_id", index = "shannon")
compareDiversity <- function(whe, groupBy, index = "shannon") {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    cd <- colData(whe)
    if (!groupBy %in% colnames(cd))
        stop("'", groupBy, "' not found in colData. ",
             "Available columns: ",
             paste(colnames(cd), collapse = ", "))

    div <- calcDiversity(whe, indices = index)
    groups <- as.character(cd[[groupBy]])
    values <- div[[index]]
    unique_groups <- unique(groups[!is.na(groups)])

    if (length(unique_groups) < 2L) {
        message("Need at least 2 groups for comparison; ",
                "returning diversity table.")
        return(div)
    }

    if (length(unique_groups) == 2L) {
        g1 <- values[groups == unique_groups[1L]]
        g2 <- values[groups == unique_groups[2L]]
        if (length(g1) < 2L || length(g2) < 2L) {
            return(DataFrame(test = "wilcoxon", statistic = NA_real_,
                             p_value = NA_real_,
                             note = "Insufficient samples per group"))
        }
        tt <- wilcox.test(g1, g2)
        DataFrame(test = "wilcoxon", statistic = tt$statistic,
                  p_value = tt$p.value)
    } else {
        tt <- kruskal.test(values ~ factor(groups))
        DataFrame(test = "kruskal", statistic = tt$statistic,
                  p_value = tt$p.value)
    }
}
