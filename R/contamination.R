# R/contamination.R
# Cross-contamination / index-hopping detection

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors DataFrame mcols
NULL

#' Detect potential cross-contamination between samples
#'
#' Identifies iSNV sites where a high-frequency variant in one
#' sample appears at suspiciously low frequency in another sample
#' from the same sequencing run, consistent with index hopping or
#' sample cross-contamination.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with multiple
#'   samples.
#' @param runCol Character scalar (optional). Column name in
#'   \code{colData} identifying the sequencing run/batch. If
#'   \code{NULL}, all samples are compared against each other. A
#'   name that is not a column of \code{colData} is an error.
#'   Samples whose run is \code{NA} are left out of every
#'   comparison, with a warning giving their number.
#' @param sourceMinFreq Numeric scalar. Minimum frequency in the
#'   "source" sample for a variant to be considered a potential
#'   contaminant origin (default: 0.20).
#' @param sinkMaxFreq Numeric scalar. Maximum frequency in the
#'   "sink" sample for a variant to be flagged as potential
#'   contamination (default: 0.05).
#' @param sinkMinFreq Numeric scalar. Minimum frequency in the
#'   "sink" sample (default: 0.005). Variants below this are
#'   likely noise, not contamination.
#' @param usePassedOnly Logical. If \code{TRUE} (default for source),
#'   only consider QC-passed variants in the source sample.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{source_sample}{Sample where the variant is at high frequency.}
#'     \item{sink_sample}{Sample where the variant appears at
#'       suspiciously low frequency.}
#'     \item{chrom}{Chromosome/segment.}
#'     \item{position}{Genomic position.}
#'     \item{ref}{Reference allele.}
#'     \item{alt}{Alternative allele.}
#'     \item{source_freq}{Frequency in the source sample.}
#'     \item{sink_freq}{Frequency in the sink sample.}
#'     \item{freq_ratio}{sink_freq / source_freq -- lower values
#'       are more suspicious.}
#'     \item{run}{Sequencing run (if \code{runCol} provided).}
#'   }
#'
#' @details
#' Index hopping occurs when adapters are exchanged between
#' library molecules during pooled sequencing, causing reads from
#' one sample to be misassigned to another. The signature: a
#' variant at high frequency (e.g., 80\%) in sample A appears at
#' very low frequency (e.g., 0.5-2\%) in sample B sequenced on
#' the same flow cell, at a ratio consistent with the expected
#' hopping rate (typically 0.1-2\% for Illumina platforms).
#'
#' This function does not require external metadata beyond what is
#' in \code{colData}. For finer control, use the output
#' \code{DataFrame} to implement study-specific contamination
#' criteria.
#'
#' @references
#' Bendall EE et al. (2022). SARS-CoV-2 genomic diversity in
#' households highlights the challenges of sequence-based
#' transmission inference. \emph{mSphere} 7:e00400-22.
#' \doi{10.1128/msphere.00400-22}
#'
#' Costello M et al. (2018). Characterization and remediation of
#' sample index swaps by non-redundant dual indexing on massively
#' parallel sequencing platforms. \emph{BMC Genomics} 19:332.
#' \doi{10.1186/s12864-018-4703-0}
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
#' ## Sample A has pos 100 at 80%; sample B has same pos at 1%
#' freq <- matrix(c(0.80, 0.05, 0.30,   # sample A
#'                  0.01, 0.40, 0.00),   # sample B
#'                nrow = 3, ncol = 2)
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freq),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("A", "B"),
#'         run = c("run1", "run1")))
#' detectCrossContamination(whe, runCol = "run")
detectCrossContamination <- function(whe,
                                     runCol = NULL,
                                     sourceMinFreq = 0.20,
                                     sinkMaxFreq = 0.05,
                                     sinkMinFreq = 0.005,
                                     usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!is.numeric(sourceMinFreq) || sourceMinFreq <= 0 ||
        sourceMinFreq >= 1)
        stop("'sourceMinFreq' must be between 0 and 1.")
    if (!is.numeric(sinkMaxFreq) || sinkMaxFreq <= 0)
        stop("'sinkMaxFreq' must be positive.")
    if (sinkMaxFreq >= sourceMinFreq)
        stop("'sinkMaxFreq' must be < 'sourceMinFreq'.")

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }

    cd <- colData(whe)
    sample_ids <- as.character(cd$sample_id)
    n_samples <- ncol(freq_mat)
    n_sites <- nrow(freq_mat)

    rr <- rowRanges(whe)
    mc <- mcols(rr)
    chr_vals <- as.character(seqnames(rr))
    pos_vals <- start(rr)
    ref_vals <- if ("ref" %in% colnames(mc)) {
        as.character(mc$ref)
    } else {
        rep(NA_character_, n_sites)
    }
    alt_vals <- if ("alt" %in% colnames(mc)) {
        as.character(mc$alt)
    } else {
        rep(NA_character_, n_sites)
    }

    ## Determine run groups
    if (is.null(runCol)) {
        runs <- rep("all", n_samples)
    } else {
        if (!is.character(runCol) || length(runCol) != 1L ||
            is.na(runCol))
            stop("'runCol' must be NULL or a single column name of ",
                 "colData(whe).")
        if (!runCol %in% colnames(cd))
            stop("Column '", runCol, "' not found in colData. ",
                 "Available columns: ",
                 paste(colnames(cd), collapse = ", "))
        runs <- as.character(cd[[runCol]])
        if (anyNA(runs))
            warning(sum(is.na(runs)), " sample(s) with NA in '", runCol,
                    "' are excluded from the comparison.")
    }

    results <- list()

    for (run_id in unique(runs)) {
        run_idx <- which(runs == run_id)
        if (length(run_idx) < 2L) next

        ## For each pair of samples in this run
        for (si in seq_along(run_idx)) {
            for (sj in seq_along(run_idx)) {
                if (si == sj) next
                src <- run_idx[si]
                snk <- run_idx[sj]

                src_freq <- freq_mat[, src]
                snk_freq <- freq_mat[, snk]

                ## Apply QC to source
                if (!is.null(qc_mat)) {
                    src_freq[!qc_mat[, src]] <- NA
                }

                ## Find sites where source is high, sink is suspicious
                suspect <- which(
                    !is.na(src_freq) & !is.na(snk_freq) &
                    src_freq >= sourceMinFreq &
                    snk_freq >= sinkMinFreq &
                    snk_freq <= sinkMaxFreq
                )

                if (length(suspect) > 0L) {
                    results <- c(results, list(data.frame(
                        source_sample = sample_ids[src],
                        sink_sample   = sample_ids[snk],
                        chrom         = chr_vals[suspect],
                        position      = pos_vals[suspect],
                        ref           = ref_vals[suspect],
                        alt           = alt_vals[suspect],
                        source_freq   = src_freq[suspect],
                        sink_freq     = snk_freq[suspect],
                        freq_ratio    = snk_freq[suspect] / src_freq[suspect],
                        run           = run_id,
                        stringsAsFactors = FALSE
                    )))
                }
            }
        }
    }

    if (length(results) == 0L) {
        return(DataFrame(
            source_sample = character(),
            sink_sample   = character(),
            chrom         = character(),
            position      = integer(),
            ref           = character(),
            alt           = character(),
            source_freq   = numeric(),
            sink_freq     = numeric(),
            freq_ratio    = numeric(),
            run           = character()
        ))
    }

    DataFrame(do.call(rbind, results))
}
