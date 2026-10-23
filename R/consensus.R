# R/consensus.R
# Multi-caller consensus for iSNV detection

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges GRanges start findOverlaps
#' @importFrom GenomeInfoDb seqnames seqinfo
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
#'   consensus sites, with a column for every sample of the inputs.
#'   \code{mcols(rowRanges)} includes:
#'   \describe{
#'     \item{n_callers}{Number of callers detecting the site; a caller
#'       that lists the same site and allele twice counts once.}
#'     \item{caller_names}{Comma-separated caller names.}
#'     \item{freq_spread}{Range of allele frequencies across callers.}
#'   }
#'   Where a site passes in more than one sample, these three columns
#'   describe one of those samples, not all of them.
#'
#'   The \code{altFreq} assay contains the mean frequency across
#'   callers for each consensus site and \code{NA} where the site did
#'   not pass in that sample, so a sample with no consensus site is a
#'   column of \code{NA}. The constructor adds a \code{qcPass} assay that
#'   is \code{TRUE} throughout. Read counts, depths and the QC history
#'   are not carried over, so run \code{\link{flagISNV}} on the
#'   individual callers' objects rather than on the consensus.
#'
#'   The \code{seqinfo()} of the result merges those of the inputs; a
#'   sequence whose length differs between inputs is an error. Each
#'   sample's \code{colData} row comes from the first input that holds
#'   the sample, so host, role and transmission-pair columns carry over.
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

    ## One row per non-missing call, ordered by caller, then sample
    ## column, then site. A record is one site and allele in one sample.
    calls <- do.call(rbind, lapply(seq_along(whe_list), function(ci)
        .consensusCalls(whe_list[[ci]], ci)))
    if (nrow(calls) == 0L) {
        stop("No variant sites found across callers.")
    }
    rec <- .firstSeenGroup(calls$chrom, calls$pos, calls$alt,
                           calls$sample_id)
    n_rec <- max(rec)

    ## Callers per record, each counted once, and the spread of their
    ## frequencies; ordering by frequency puts each record's minimum
    ## first and its maximum last
    once <- !duplicated((rec - 1) * length(whe_list) + calls$caller)
    n_callers <- tabulate(rec[once], n_rec)
    o <- order(rec, calls$freq)
    spread <- calls$freq[o][!duplicated(rec[o], fromLast = TRUE)] -
        calls$freq[o][!duplicated(rec[o])]

    kept <- which(n_callers >= min_callers & spread <= freq_tolerance)
    if (length(kept) == 0L) {
        stop("No sites passed consensus filters ",
             "(min_callers = ", min_callers,
             ", freq_tolerance = ", freq_tolerance, ").")
    }
    in_kept <- rec %in% kept
    mean_freq <- vapply(
        split(calls$freq[in_kept], factor(rec[in_kept], levels = kept)),
        mean, numeric(1), USE.NAMES = FALSE)
    names_kept <- .callerNames(rec[once], calls$caller[once], kept,
                               caller_names)

    ## Each kept record by its first call; sites in the order they are
    ## first seen among the kept records, and every sample of the inputs
    first <- which(!duplicated(rec))[kept]
    site <- .firstSeenGroup(calls$chrom[first], calls$pos[first],
                            calls$alt[first])
    input_ids <- lapply(whe_list, .consensusSampleIds)
    all_samples <- unique(unlist(input_ids, use.names = FALSE))
    freq_mat <- matrix(NA_real_, nrow = max(site),
                       ncol = length(all_samples),
                       dimnames = list(NULL, all_samples))
    freq_mat[cbind(site, match(calls$sample_id[first], all_samples))] <-
        mean_freq

    ## Coordinates come from each site's own calls. The site columns
    ## describe the last kept record of each site, as they always have.
    at <- first[!duplicated(site)]
    last <- integer(max(site))
    last[site] <- seq_along(site)
    gr <- GRanges(seqnames = calls$chrom[at],
                  ranges = IRanges(calls$pos[at], width = 1L),
                  seqinfo = .mergeSeqinfo(whe_list))
    mcols(gr)$ref <- calls$ref[first][last]
    mcols(gr)$alt <- calls$alt[first][last]
    mcols(gr)$n_callers <- n_callers[kept][last]
    mcols(gr)$caller_names <- names_kept[last]
    mcols(gr)$freq_spread <- spread[kept][last]

    cd <- .consensusColData(whe_list, input_ids, all_samples)

    WithinHostExperiment(
        assays = list(altFreq = freq_mat),
        rowRanges = gr,
        colData = cd
    )
}

## The non-missing calls of one caller's object, one row each, sample
## column by sample column; 'ci' is the caller's place in the list.
.consensusCalls <- function(w, ci) {
    rr <- rowRanges(w)
    mc <- mcols(rr)
    fm <- as.matrix(assay(w, "altFreq"))
    sample_ids <- .consensusSampleIds(w)
    allele <- function(col) {
        if (col %in% colnames(mc)) as.character(mc[[col]])
        else rep(NA_character_, length(rr))
    }
    hit <- which(!is.na(fm), arr.ind = TRUE)
    i <- hit[, 1L]
    data.frame(chrom = as.character(seqnames(rr))[i],
               pos = start(rr)[i],
               ref = allele("ref")[i],
               alt = allele("alt")[i],
               sample_id = sample_ids[hit[, 2L]],
               caller = rep(ci, length(i)),
               freq = fm[hit],
               stringsAsFactors = FALSE)
}

## Sample identifiers of one caller's object: its column names, or S1,
## S2, ... when it has none.
.consensusSampleIds <- function(w) {
    ids <- colnames(w)
    if (is.null(ids)) paste0("S", seq_len(ncol(w))) else ids
}

## The inputs' seqinfo merged into one. A sequence whose length or
## genome differs between inputs is an error; inputs that list
## different sequences are combined without the merge's warning.
.mergeSeqinfo <- function(whe_list) {
    infos <- lapply(whe_list, seqinfo)
    tryCatch(
        withCallingHandlers(
            Reduce(merge, infos),
            warning = function(w) {
                if (grepl("sequence levels", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }),
        error = function(e) {
            stop("The objects in 'whe_list' disagree on their sequences: ",
                 conditionMessage(e), call. = FALSE)
        })
}

## colData for the consensus: each sample's row from the first input
## that holds the sample, with the columns of those inputs combined and
## NA where an input lacks a column. 'ids' lists each input's samples.
.consensusColData <- function(whe_list, ids, samples) {
    src <- integer(length(samples))
    row <- integer(length(samples))
    for (i in rev(seq_along(whe_list))) {
        m <- match(samples, ids[[i]])
        src[!is.na(m)] <- i
        row[!is.na(m)] <- m[!is.na(m)]
    }
    used <- unique(src)
    cds <- lapply(whe_list[used], colData)
    cols <- unique(unlist(lapply(cds, colnames), use.names = FALSE))
    pieces <- lapply(seq_along(used), function(k) {
        piece <- cds[[k]][row[src == used[k]], , drop = FALSE]
        for (col in setdiff(cols, colnames(piece))) piece[[col]] <- NA
        piece[, cols, drop = FALSE]
    })
    out <- do.call(rbind, unname(pieces))
    placed <- unlist(lapply(used, function(i) which(src == i)))
    out <- out[order(placed), , drop = FALSE]
    out$sample_id <- samples
    rownames(out) <- samples
    out
}

## Number the distinct combinations of the vectors in '...' (all of one
## length) in the order they are first seen, with no pasted key, so no
## character inside a value can make two combinations collide.
.firstSeenGroup <- function(...) {
    id <- rep.int(1L, length(..1))
    for (x in list(...)) {
        k <- match(x, unique(x))
        key <- (as.numeric(id) - 1) * max(k) + k
        id <- match(key, unique(key))
    }
    id
}

## Comma-separated caller names of each kept record, in list order.
## 'rec' and 'caller' hold one entry per caller of a record.
.callerNames <- function(rec, caller, kept, caller_names) {
    at <- match(rec, kept)
    out <- character(length(kept))
    started <- logical(length(kept))
    for (ci in seq_along(caller_names)) {
        j <- at[caller == ci & !is.na(at)]
        out[j] <- ifelse(started[j], paste0(out[j], ",", caller_names[ci]),
                         caller_names[ci])
        started[j] <- TRUE
    }
    out
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
