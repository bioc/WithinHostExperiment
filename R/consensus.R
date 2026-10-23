# R/consensus.R
# Multi-caller consensus for iSNV detection

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges GRanges start findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom methods is new
NULL

# ============================================================
# consensusISNV -- multi-caller consensus
# ============================================================

#' Multi-caller consensus iSNV detection
#'
#' Takes a named list of \code{\link{WithinHostExperiment}} objects
#' (one per variant caller) and identifies variant sites detected
#' by multiple callers. Sites are matched by genomic position
#' (seqnames + start) and alternative allele. Only sites detected
#' by at least \code{min_callers} with a frequency spread no
#' greater than \code{freq_tolerance} are retained.
#'
#' When technical replicates are unavailable, using a multi-caller
#' consensus is recommended for reducing false positive iSNV calls
#' (Roder et al. 2023).
#'
#' @param whe_list A named list of \code{\link{WithinHostExperiment}}
#'   objects, one per variant caller (e.g.,
#'   \code{list(ivar = whe1, lofreq = whe2)}).
#' @param min_callers Integer scalar. Minimum number of callers
#'   that must detect a site (default: 2L).
#' @param freq_tolerance Numeric scalar. Maximum allowed range
#'   (max - min) of allele frequencies across callers for a site
#'   to pass (default: 0.05).
#'
#' @return A \code{\link{WithinHostExperiment}} containing only
#'   consensus sites. \code{mcols(rowRanges)} includes:
#'   \describe{
#'     \item{n_callers}{Number of callers detecting the site; a caller
#'       that lists the same site and allele twice counts once.}
#'     \item{caller_names}{Comma-separated caller names.}
#'     \item{freq_spread}{Range of allele frequencies across callers.}
#'   }
#'   The \code{altFreq} assay contains the mean frequency across
#'   callers for each consensus site. It is the only assay of the
#'   result: read counts, depths and \code{qcPass} are not carried
#'   over, so run \code{\link{flagISNV}} on the individual callers'
#'   objects rather than on the consensus.
#'
#' @references
#' Roder AE et al. (2023). Optimized quantification of intra-host
#' viral diversity in SARS-CoV-2 and influenza virus sequence data.
#' \emph{mBio} 14:e01046-23.
#' \doi{10.1128/mbio.01046-23}
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
#' whe1 <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.10, 0.20, 0.30), ncol = 1)),
#'     rowRanges = gr, colData = S4Vectors::DataFrame(sample_id = "S1"))
#' gr2 <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 400), width = 1))
#' S4Vectors::mcols(gr2)$ref <- c("A", "C", "T")
#' S4Vectors::mcols(gr2)$alt <- c("T", "G", "C")
#' whe2 <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.12, 0.18, 0.15), ncol = 1)),
#'     rowRanges = gr2, colData = S4Vectors::DataFrame(sample_id = "S1"))
#' consensus <- consensusISNV(list(caller1 = whe1, caller2 = whe2))
#' consensus
consensusISNV <- function(whe_list, min_callers = 2L,
                          freq_tolerance = 0.05) {
    if (!is.list(whe_list) || length(whe_list) < 2L)
        stop("'whe_list' must be a list of at least 2 ",
             "WithinHostExperiment objects.")
    if (is.null(names(whe_list)))
        stop("'whe_list' must be a named list (names = caller names).")
    are_whe <- vapply(whe_list, is, logical(1), "WithinHostExperiment")
    if (!all(are_whe))
        stop("All elements of 'whe_list' must be ",
             "WithinHostExperiment objects.")
    if (!is.numeric(min_callers) || length(min_callers) != 1L ||
        min_callers < 1L)
        stop("'min_callers' must be a single positive integer.")
    if (!is.numeric(freq_tolerance) || length(freq_tolerance) != 1L ||
        freq_tolerance < 0)
        stop("'freq_tolerance' must be a single non-negative number.")

    caller_names <- names(whe_list)

    ## Build a unified site table: chrom:pos:alt -> list of callers + freqs
    ## Use the first sample column from each WHE for frequency comparison
    site_records <- list()
    for (ci in seq_along(whe_list)) {
        w <- whe_list[[ci]]
        rr <- rowRanges(w)
        freq_mat <- assay(w, "altFreq")
        mc <- mcols(rr)
        chr_vals <- as.character(seqnames(rr))
        pos_vals <- start(rr)
        alt_vals <- if ("alt" %in% colnames(mc)) {
            as.character(mc$alt)
        } else {
            rep(NA_character_, length(rr))
        }
        ref_vals <- if ("ref" %in% colnames(mc)) {
            as.character(mc$ref)
        } else {
            rep(NA_character_, length(rr))
        }

        n_samples <- ncol(freq_mat)
        for (si in seq_len(n_samples)) {
            sample_id <- colnames(freq_mat)[si]
            if (is.null(sample_id)) sample_id <- paste0("S", si)
            for (ri in seq_len(length(rr))) {
                freq_val <- freq_mat[ri, si]
                if (is.na(freq_val)) next
                key <- paste(chr_vals[ri], pos_vals[ri],
                             alt_vals[ri], sample_id, sep = ":")
                if (is.null(site_records[[key]])) {
                    site_records[[key]] <- list(
                        chrom = chr_vals[ri],
                        pos = pos_vals[ri],
                        ref = ref_vals[ri],
                        alt = alt_vals[ri],
                        sample_id = sample_id,
                        callers = character(),
                        freqs = numeric()
                    )
                }
                site_records[[key]]$callers <- c(
                    site_records[[key]]$callers,
                    caller_names[ci])
                site_records[[key]]$freqs <- c(
                    site_records[[key]]$freqs, freq_val)
            }
        }
    }

    if (length(site_records) == 0L) {
        stop("No variant sites found across callers.")
    }

    ## Filter by min_callers and freq_tolerance
    keep <- vapply(site_records, function(rec) {
        n <- length(unique(rec$callers))
        if (n < min_callers) return(FALSE)
        spread <- max(rec$freqs) - min(rec$freqs)
        spread <= freq_tolerance
    }, logical(1))

    kept <- site_records[keep]
    if (length(kept) == 0L) {
        stop("No sites passed consensus filters ",
             "(min_callers = ", min_callers,
             ", freq_tolerance = ", freq_tolerance, ").")
    }

    ## Group by sample_id and build matrices
    all_samples <- unique(vapply(kept, `[[`, character(1), "sample_id"))
    ## Get unique site keys (chrom:pos:alt)
    site_keys <- unique(vapply(kept, function(rec) {
        paste(rec$chrom, rec$pos, rec$alt, sep = ":")
    }, character(1)))

    ## Build GRanges for unique sites
    site_parts <- strsplit(site_keys, ":")
    gr <- GRanges(
        seqnames = vapply(site_parts, `[`, character(1), 1L),
        ranges = IRanges(
            as.integer(vapply(site_parts, `[`, character(1), 2L)),
            width = 1L)
    )

    ## Get ref/alt/metadata from first record of each site
    n_callers_vec <- integer(length(site_keys))
    caller_names_vec <- character(length(site_keys))
    freq_spread_vec <- numeric(length(site_keys))
    ref_vec <- character(length(site_keys))
    alt_vec <- character(length(site_keys))

    n_sites <- length(site_keys)
    n_samp <- length(all_samples)
    freq_mat <- matrix(NA_real_, nrow = n_sites, ncol = n_samp,
                       dimnames = list(NULL, all_samples))

    for (rec in kept) {
        sk <- paste(rec$chrom, rec$pos, rec$alt, sep = ":")
        si <- match(sk, site_keys)
        sj <- match(rec$sample_id, all_samples)
        freq_mat[si, sj] <- mean(rec$freqs)
        n_callers_vec[si] <- length(unique(rec$callers))
        caller_names_vec[si] <- paste(unique(rec$callers), collapse = ",")
        freq_spread_vec[si] <- max(rec$freqs) - min(rec$freqs)
        ref_vec[si] <- rec$ref
        alt_vec[si] <- rec$alt
    }

    mcols(gr)$ref <- ref_vec
    mcols(gr)$alt <- alt_vec
    mcols(gr)$n_callers <- n_callers_vec
    mcols(gr)$caller_names <- caller_names_vec
    mcols(gr)$freq_spread <- freq_spread_vec

    cd <- DataFrame(sample_id = all_samples)

    WithinHostExperiment(
        assays = list(altFreq = freq_mat),
        rowRanges = gr,
        colData = cd
    )
}


# ============================================================
# calcCallerConcordance -- pairwise caller agreement
# ============================================================

#' Pairwise caller concordance matrix
#'
#' For a list of \code{\link{WithinHostExperiment}} objects (one per
#' variant caller), computes pairwise concordance statistics
#' including shared sites, caller-unique sites, and Jaccard index.
#'
#' @param whe_list A named list of \code{\link{WithinHostExperiment}}
#'   objects, one per variant caller.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{caller1}{First caller name.}
#'     \item{caller2}{Second caller name.}
#'     \item{n_shared}{Number of sites detected by both callers.}
#'     \item{n_only1}{Sites unique to caller1.}
#'     \item{n_only2}{Sites unique to caller2.}
#'     \item{jaccard}{Jaccard similarity index (shared / union).}
#'   }
#'
#' @references
#' Roder AE et al. (2023). Optimized quantification of intra-host
#' viral diversity in SARS-CoV-2 and influenza virus sequence data.
#' \emph{mBio} 14:e01046-23.
#' \doi{10.1128/mbio.01046-23}
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
#' gr2 <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 400), width = 1))
#' S4Vectors::mcols(gr2)$ref <- c("A", "C", "T")
#' S4Vectors::mcols(gr2)$alt <- c("T", "G", "C")
#' whe2 <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.12, 0.22, 0.15), ncol = 1)),
#'     rowRanges = gr2,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' calcCallerConcordance(list(ivar = whe1, lofreq = whe2))
calcCallerConcordance <- function(whe_list) {
    if (!is.list(whe_list) || length(whe_list) < 2L)
        stop("'whe_list' must be a list of at least 2 ",
             "WithinHostExperiment objects.")
    if (is.null(names(whe_list)))
        stop("'whe_list' must be a named list (names = caller names).")
    are_whe <- vapply(whe_list, is, logical(1), "WithinHostExperiment")
    if (!all(are_whe))
        stop("All elements of 'whe_list' must be ",
             "WithinHostExperiment objects.")

    caller_names <- names(whe_list)

    ## Build site key sets per caller
    site_sets <- lapply(whe_list, function(w) {
        rr <- rowRanges(w)
        mc <- mcols(rr)
        chr_vals <- as.character(seqnames(rr))
        pos_vals <- start(rr)
        alt_vals <- if ("alt" %in% colnames(mc)) {
            as.character(mc$alt)
        } else {
            rep(NA_character_, length(rr))
        }
        ## Only include sites with non-NA frequency in at least one sample
        freq_mat <- assay(w, "altFreq")
        has_data <- rowSums(!is.na(freq_mat)) > 0L
        paste(chr_vals[has_data], pos_vals[has_data],
              alt_vals[has_data], sep = ":")
    })

    ## Pairwise comparisons
    n_callers <- length(caller_names)
    pairs <- expand.grid(i = seq_len(n_callers),
                         j = seq_len(n_callers),
                         stringsAsFactors = FALSE)
    pairs <- pairs[pairs$i < pairs$j, , drop = FALSE]

    results <- lapply(seq_len(nrow(pairs)), function(p) {
        i <- pairs$i[p]
        j <- pairs$j[p]
        set_i <- site_sets[[i]]
        set_j <- site_sets[[j]]
        shared <- length(intersect(set_i, set_j))
        only_i <- length(setdiff(set_i, set_j))
        only_j <- length(setdiff(set_j, set_i))
        union_n <- shared + only_i + only_j
        jaccard <- if (union_n > 0L) shared / union_n else NA_real_
        data.frame(
            caller1 = caller_names[i],
            caller2 = caller_names[j],
            n_shared = shared,
            n_only1 = only_i,
            n_only2 = only_j,
            jaccard = jaccard,
            stringsAsFactors = FALSE)
    })

    DataFrame(do.call(rbind, results))
}
