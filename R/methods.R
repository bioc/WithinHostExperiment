# R/methods.R
# Constructor functions, accessor methods, show methods, internal helpers

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods setMethod setAs new callNextMethod is show
#' @importFrom BiocGenerics nrow ncol
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom VariantAnnotation VRanges
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   colData rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
NULL

# ============================================================
# ConcordanceRule constructor + show
# ============================================================

#' Create a ConcordanceRule object
#'
#' Says how the technical replicates of a sample must agree, for use as
#' the \code{replicateConc} of an \code{\link{ISNVFilter}}.
#'
#' @param freqTolerance Numeric. The largest difference in allele
#'   frequency allowed between two replicates that both made the call
#'   (default: 0.05); \code{Inf} compares no frequencies.
#' @param requireBothDetected Logical. Whether a site called in one
#'   replicate of a pair but not in the other fails (default:
#'   \code{TRUE}).
#'
#' @details Neither argument is converted from another type: a factor, a
#'   string or a logical \code{freqTolerance} is an error, and so is a
#'   \code{requireBothDetected} that is not \code{TRUE} or \code{FALSE}, as
#'   in \code{\link{flagReplicateDiscordance}}. An integer tolerance is
#'   stored as a double.
#'
#' @return A \code{\linkS4class{ConcordanceRule}} object.
#' @seealso \code{\link{flagISNV}}, which applies the rule, and
#'   \code{\link{flagReplicateDiscordance}}.
#' @export
#' @examples
#' rule <- ConcordanceRule(freqTolerance = 0.05)
#' rule
#' ISNVFilter(replicateConc = rule)
ConcordanceRule <- function(freqTolerance = 0.05,
                            requireBothDetected = TRUE) {
    ## the slots would store a factor's level codes, so check the type
    ## here rather than convert it
    if (!is.numeric(freqTolerance) || is.factor(freqTolerance))
        stop("'freqTolerance' must be a single positive number, not ",
             class(freqTolerance)[1L], ".")
    if (!is.logical(requireBothDetected))
        stop("'requireBothDetected' must be TRUE or FALSE, not ",
             class(requireBothDetected)[1L], ".")
    new("ConcordanceRule",
        freqTolerance       = as.numeric(freqTolerance),
        requireBothDetected = requireBothDetected)
}

#' @rdname ConcordanceRule-class
#' @param object A \code{ConcordanceRule} object.
#' @export
setMethod("show", "ConcordanceRule", function(object) {
    cat("ConcordanceRule object\n")
    cat("  freqTolerance:      ", slot(object, "freqTolerance"), "\n")
    cat("  requireBothDetected:", slot(object, "requireBothDetected"),
        "\n")
})


# ============================================================
# ISNVFilter constructor + show
# ============================================================

#' Create an ISNVFilter object
#'
#' Constructs a quality control filter configuration for use with
#' \code{\link{flagISNV}}.
#'
#' @param minDepth Integer. Minimum read depth (default: 100).
#' @param minFreq Numeric. Minimum alt allele frequency (default: 0.03).
#' @param maxFreq Numeric. Maximum alt allele frequency (default: 0.50).
#' @param minAltReads Integer. Minimum alt allele reads (default: 10).
#' @param maxStrandBias Numeric. Max strand bias fold-diff (default: 10).
#' @param replicateConc A \code{\link{ConcordanceRule}} saying how
#'   technical replicates must agree, which \code{\link{flagISNV}}
#'   applies after the other thresholds, or \code{NULL} (default) for
#'   none.
#'
#' @return An \code{\link{ISNVFilter}} object.
#' @export
#' @examples
#' ISNVFilter()
#' ISNVFilter(minDepth = 500L, minFreq = 0.05)
#' ISNVFilter(replicateConc = ConcordanceRule(freqTolerance = 0.05))
ISNVFilter <- function(minDepth = 100L, minFreq = 0.03, maxFreq = 0.50,
                       minAltReads = 10L, maxStrandBias = 10.0,
                       replicateConc = NULL) {
    if (!is.null(replicateConc) && !is(replicateConc, "ConcordanceRule"))
        stop("'replicateConc' must be a ConcordanceRule or NULL; ",
             "replicateConc = ConcordanceRule() requires technical ",
             "replicates to agree.")
    new("ISNVFilter",
        minDepth      = as.integer(minDepth),
        minFreq       = as.numeric(minFreq),
        maxFreq       = as.numeric(maxFreq),
        minAltReads   = as.integer(minAltReads),
        maxStrandBias = as.numeric(maxStrandBias),
        replicateConc = replicateConc)
}

#' @rdname ISNVFilter-class
#' @param object An \code{ISNVFilter} object.
#' @export
setMethod("show", "ISNVFilter", function(object) {
    cat("ISNVFilter object\n")
    cat("  minDepth:     ", slot(object, "minDepth"), "\n")
    cat("  freq range:   [", slot(object, "minFreq"), ", ",
        slot(object, "maxFreq"), "]\n", sep = "")
    cat("  minAltReads:  ", slot(object, "minAltReads"), "\n")
    cat("  maxStrandBias:", slot(object, "maxStrandBias"), "\n")
    rule <- slot(object, "replicateConc")
    if (is.null(rule)) {
        cat("  replicateConc: none\n")
    } else {
        cat("  replicateConc: freqTolerance ",
            slot(rule, "freqTolerance"), ", requireBothDetected ",
            slot(rule, "requireBothDetected"), "\n", sep = "")
    }
})


# ============================================================
# WithinHostExperiment constructor
# ============================================================

#' Create a WithinHostExperiment from components
#'
#' Constructs a \code{\link{WithinHostExperiment}} from assay matrices,
#' genomic ranges, and sample metadata.
#'
#' @param assays Named list of matrices. Must include at least
#'   \code{altFreq}. If \code{qcPass} is not provided, it is
#'   automatically created as an all-TRUE matrix.
#' @param rowRanges \code{GRanges}. Genomic coordinates for each
#'   variant site.
#' @param colData \code{DataFrame} or \code{data.frame}. Sample
#'   metadata. Must contain a \code{sample_id} column.
#' @param metadata Named list (optional). Additional metadata such as
#'   \code{reference_genome}, \code{genome_length}, \code{calling_info}.
#'
#' @return A \code{\link{WithinHostExperiment}} object.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("seg1", IRanges::IRanges(c(100, 200, 300), width = 1))
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' whe
WithinHostExperiment <- function(assays, rowRanges, colData,
                                 metadata = list()) {
    ## Coerce colData
    if (is.data.frame(colData) && !is(colData, "DataFrame")) {
        colData <- DataFrame(colData)
    }

    ## Validate sample_id
    if (!"sample_id" %in% colnames(colData)) {
        stop("'colData' must contain a 'sample_id' column. ",
             "Found columns: ",
             paste(colnames(colData), collapse = ", "))
    }

    ## Every downstream function reads altFreq
    if (!is.list(assays) || !length(assays) ||
            !"altFreq" %in% names(assays)) {
        stop("'assays' must be a named list containing 'altFreq'. ",
             "Found: ",
             if (length(names(assays))) {
                 paste(names(assays), collapse = ", ")
             } else {
                 "no named assays"
             })
    }

    n_sites   <- nrow(assays[[1]])
    n_samples <- ncol(assays[[1]])
    sample_ids <- as.character(colData$sample_id)

    ## Set column names on assay matrices from sample_id
    for (nm in names(assays)) {
        if (is.null(colnames(assays[[nm]]))) {
            colnames(assays[[nm]]) <- sample_ids
        }
    }

    ## Auto-create qcPass if missing
    if (!"qcPass" %in% names(assays)) {
        assays[["qcPass"]] <- matrix(TRUE,
            nrow = n_sites, ncol = n_samples,
            dimnames = list(NULL, sample_ids))
    }

    ## Build SE
    se <- SummarizedExperiment(
        assays    = assays,
        rowRanges = rowRanges,
        colData   = colData
    )
    S4Vectors::metadata(se) <- metadata

    ## Construct WHE
    new("WithinHostExperiment",
        se,
        qcLog = DataFrame(
            step      = character(),
            parameter = character(),
            value     = character(),
            n_flagged = integer(),
            n_passed  = integer(),
            timestamp = character()
        ),
        qcFilters = list()
    )
}


# ============================================================
# WithinHostExperiment accessor methods
# ============================================================

#' @rdname WHE-accessors
#' @export
setMethod("qcLog", "WithinHostExperiment", function(x) x@qcLog)

#' @rdname WHE-accessors
#' @export
setMethod("qcFilters", "WithinHostExperiment", function(x) x@qcFilters)

#' @rdname WHE-accessors
#' @export
setReplaceMethod("qcLog", "WithinHostExperiment", function(x, value) {
    x@qcLog <- value
    x
})

#' @rdname WHE-accessors
#' @export
setReplaceMethod("qcFilters", "WithinHostExperiment", function(x, value) {
    x@qcFilters <- value
    x
})

#' @rdname WHE-accessors
#' @export
setMethod("passedISNV", "WithinHostExperiment", function(x) {
    if (!"qcPass" %in% assayNames(x)) {
        return(x)
    }
    qc_mat <- assay(x, "qcPass")
    ## Keep rows where at least one sample passed
    keep <- rowSums(qc_mat, na.rm = TRUE) > 0L
    x[keep, ]
})

#' @rdname WHE-accessors
#' @export
setMethod("transmissionPairs", "WithinHostExperiment", function(x) {
    cd <- colData(x)

    ## Check required columns exist
    if (!"role" %in% colnames(cd) || !"pair_id" %in% colnames(cd)) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }

    ## Filter to donor/recipient rows with non-NA pair_id
    has_role <- !is.na(cd$role) & cd$role %in% c("donor", "recipient")
    has_pair <- !is.na(cd$pair_id)
    relevant <- cd[has_role & has_pair, , drop = FALSE]

    if (nrow(relevant) == 0L) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }

    ## Group by pair_id
    pair_ids <- unique(as.character(relevant$pair_id))
    result <- lapply(pair_ids, function(pid) {
        sub       <- relevant[relevant$pair_id == pid, , drop = FALSE]
        donors    <- as.character(sub$sample_id[sub$role == "donor"])
        recipients <- as.character(sub$sample_id[sub$role == "recipient"])
        if (length(donors) == 0L || length(recipients) == 0L) {
            return(NULL)
        }
        expand.grid(donor     = donors,
                    recipient = recipients,
                    pair_id   = pid,
                    stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, Filter(Negate(is.null), result))
    if (is.null(result) || nrow(result) == 0L) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }
    DataFrame(result)
})


# ============================================================
# WithinHostExperiment show method
# ============================================================

#' @rdname WithinHostExperiment-class
#' @param object A \code{WithinHostExperiment} object.
#' @export
setMethod("show", "WithinHostExperiment", function(object) {
    cat("WithinHostExperiment with",
        nrow(object), "variant sites and",
        ncol(object), "samples\n")
    anames <- assayNames(object)
    if (length(anames) > 0L) {
        cat("  assays(", length(anames), "): ",
            paste(anames, collapse = ", "), "\n", sep = "")
    }
    cat("  qcLog:", nrow(qcLog(object)), "filtering step(s)\n")
    pairs <- transmissionPairs(object)
    cat("  transmission pairs:", nrow(pairs), "\n")
})


# ============================================================
# Internal helper functions
# ============================================================

#' Record a QC step in the audit log
#'
#' @param whe A \code{WithinHostExperiment} object.
#' @param step Character. Filter step name.
#' @param parameter Character. Parameter name.
#' @param value Character. Parameter value (as string).
#' @param n_flagged Integer. Number of entries newly flagged.
#' @param n_passed Integer. Number of entries still passing.
#' @return Updated \code{WithinHostExperiment}.
#' @keywords internal
#' @noRd
.addQcStep <- function(whe, step, parameter, value, n_flagged, n_passed) {
    entry <- DataFrame(
        step      = step,
        parameter = parameter,
        value     = as.character(value),
        n_flagged = as.integer(n_flagged),
        n_passed  = as.integer(n_passed),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )
    qcLog(whe) <- rbind(qcLog(whe), entry)
    whe
}

## TRUE where a qcPass cell passed; an NA cell counts as not passed. Use
## as freq[!.passedQC(qc)] <- NA: freq[!qc] <- NA skips the NA cells.
.passedQC <- function(qc) !is.na(qc) & qc

#' Internal constructor from pre-validated components
#'
#' @param rowRanges \code{GRanges}. Variant positions.
#' @param assays Named list of matrices.
#' @param colData \code{DataFrame}. Sample metadata.
#' @param callerName Character. Caller name for metadata.
#' @param metadata Named list. Additional metadata.
#' @return A \code{WithinHostExperiment}.
#' @keywords internal
#' @noRd
.makeWHE <- function(rowRanges, assays, colData,
                     callerName = NA_character_, metadata = list()) {
    sample_ids <- as.character(colData$sample_id)
    n_sites    <- length(rowRanges)
    n_samples  <- length(sample_ids)

    ## Ensure colnames on all assays
    for (nm in names(assays)) {
        if (is.null(colnames(assays[[nm]]))) {
            colnames(assays[[nm]]) <- sample_ids
        }
    }

    ## Auto-create qcPass
    if (!"qcPass" %in% names(assays)) {
        assays[["qcPass"]] <- matrix(TRUE,
            nrow = n_sites, ncol = n_samples,
            dimnames = list(NULL, sample_ids))
    }

    ## Build SE
    se <- SummarizedExperiment(
        assays    = assays,
        rowRanges = rowRanges,
        colData   = colData
    )

    ## Store caller in metadata
    if (!is.na(callerName)) {
        metadata[["calling_info"]] <- DataFrame(
            caller  = callerName,
            version = NA_character_
        )
    }
    S4Vectors::metadata(se) <- metadata

    new("WithinHostExperiment",
        se,
        qcLog = DataFrame(
            step      = character(),
            parameter = character(),
            value     = character(),
            n_flagged = integer(),
            n_passed  = integer(),
            timestamp = character()
        ),
        qcFilters = list()
    )
}


# ============================================================
# Subsetting -- preserve custom slots
# ============================================================

#' Subset a WithinHostExperiment
#'
#' Subsetting preserves the \code{qcLog} and \code{qcFilters} slots.
#' These slots are experiment-level metadata (not row/column-dependent)
#' and are carried forward unchanged.
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param i Row indices (variant sites).
#' @param j Column indices (samples).
#' @param ... Additional arguments passed to the parent method.
#' @param drop Logical (ignored for \code{SummarizedExperiment}).
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @rdname WithinHostExperiment-class
#' @export
setMethod("[", "WithinHostExperiment",
    function(x, i, j, ..., drop = TRUE) {
        result <- callNextMethod()
        new("WithinHostExperiment", result,
            qcLog = qcLog(x),
            qcFilters = qcFilters(x))
    }
)

# ============================================================
# Combine WithinHostExperiment objects
# ============================================================

#' Combine WithinHostExperiment objects by rows or columns
#'
#' Methods for the \code{\link[S4Vectors]{combineRows}} and
#' \code{\link[S4Vectors]{combineCols}} generics.
#' \code{combineRows} combines objects vertically (adding variant
#' sites); all objects must have the same samples (columns). Samples
#' are matched by their column names: an object that holds them in
#' another order is put in the first object's order, and one that holds
#' other samples is an error naming them. The \code{colData} and the
#' \code{metadata} are the first object's.
#' \code{combineCols} combines objects horizontally (adding
#' samples); all objects must have the same variant sites (rows).
#' Sites are matched by sequence name, position, \code{ref} and
#' \code{alt}: an object that holds them in another order is put in the
#' first object's order, and one that holds other sites is an error
#' counting them. Objects read from separate files usually hold
#' different sites; read the files together with
#' \code{\link{readWithinHost}} instead, which takes the union.
#' QC logs and filters from all objects are merged.
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param ... Further \code{\link{WithinHostExperiment}} objects to
#'   combine.
#' @param use.names Accepted for compatibility with the generic. Sites
#'   are matched by sequence name, position and alleles, not by row
#'   name.
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @name WHE-combine
#' @aliases combineRows,WithinHostExperiment-method
#'   combineCols,WithinHostExperiment-method
#' @importFrom S4Vectors combineRows combineCols
#' @examples
#' data(example_whe)
#' whe1 <- example_whe[1:10, ]
#' whe2 <- example_whe[11:nrow(example_whe), ]
#' merged <- combineRows(whe1, whe2)
#' nrow(merged) == nrow(example_whe)
#'
#' whe_d <- example_whe[, 1]
#' whe_r <- example_whe[, 2]
#' merged <- combineCols(whe_d, whe_r)
#' ncol(merged)
NULL

#' @rdname WHE-combine
#' @export
setMethod("combineRows", "WithinHostExperiment", function(x, ...) {
    args <- c(list(x), list(...))
    if (!all(vapply(args, is, logical(1), "WithinHostExperiment")))
        stop("All arguments must be WithinHostExperiment objects.")
    ## Samples are matched by name: the same samples in another order are
    ## put in the first object's order, and other samples are an error
    args <- .alignTo(args, .colKey, "samples", "combineRows()")
    ## Merge qcLogs and qcFilters
    all_logs <- .mergeQcLogs(args)
    all_filters <- unlist(lapply(args, qcFilters),
                          recursive = FALSE)

    ## Concatenate rowRanges
    all_rr <- do.call(c, lapply(args, rowRanges))

    ## Get common colData from first object
    cd <- colData(args[[1L]])

    ## Concatenate assay matrices by rows
    anames <- assayNames(args[[1L]])
    merged_assays <- lapply(anames, function(nm) {
        mats <- lapply(args, function(obj) assay(obj, nm))
        do.call(base::rbind, mats)
    })
    names(merged_assays) <- anames

    ## Build merged WHE
    se <- SummarizedExperiment(
        assays = merged_assays, rowRanges = all_rr, colData = cd)
    S4Vectors::metadata(se) <- S4Vectors::metadata(args[[1L]])
    new("WithinHostExperiment", se,
        qcLog = all_logs, qcFilters = all_filters)
})

#' @rdname WHE-combine
#' @export
setMethod("combineCols", "WithinHostExperiment",
    function(x, ..., use.names = TRUE) {
    args <- c(list(x), list(...))
    if (!all(vapply(args, is, logical(1), "WithinHostExperiment")))
        stop("All arguments must be WithinHostExperiment objects.")
    ## Sites are matched by sequence, position and alleles: the same sites
    ## in another order are put in the first object's order, and other
    ## sites are an error
    args <- .alignTo(args, .rowKey, "sites", "combineCols()", paste0(
        " Objects read from separate files usually hold different sites; ",
        "read the files together, as readWithinHost(c(vcf1, vcf2)) does, ",
        "to take the union of their sites."))
    all_logs <- .mergeQcLogs(args)
    all_filters <- unlist(lapply(args, qcFilters),
                          recursive = FALSE)

    ## Common rowRanges from first object -- strip all names
    ## to avoid SummarizedExperiment cbind validation conflicts
    rr <- rowRanges(args[[1L]])
    names(rr) <- NULL
    rownames(rr) <- NULL
    if (ncol(S4Vectors::mcols(rr)) > 0L) {
        rownames(S4Vectors::mcols(rr)) <- NULL
    }

    ## Concatenate colData
    all_cd <- do.call(
        function(...) base::rbind(...),
        lapply(args, function(obj) as.data.frame(colData(obj))))
    all_cd <- DataFrame(all_cd)
    rownames(all_cd) <- NULL

    ## Concatenate assay matrices by columns
    anames <- assayNames(args[[1L]])
    merged_assays <- lapply(anames, function(nm) {
        mats <- lapply(args, function(obj) {
            m <- assay(obj, nm)
            rownames(m) <- NULL
            m
        })
        do.call(base::cbind, mats)
    })
    names(merged_assays) <- anames

    se <- SummarizedExperiment(
        assays = merged_assays, rowRanges = rr, colData = all_cd)
    S4Vectors::metadata(se) <- S4Vectors::metadata(args[[1L]])
    new("WithinHostExperiment", se,
        qcLog = all_logs, qcFilters = all_filters)
})

## The seqname:start:ref:alt keys by which combineCols() matches rows
.rowKey <- function(obj) {
    rr <- rowRanges(obj)
    mc <- S4Vectors::mcols(rr)
    allele <- function(nm) if (nm %in% colnames(mc)) as.character(mc[[nm]])
                           else rep(NA_character_, length(rr))
    paste(as.character(GenomeInfoDb::seqnames(rr)), BiocGenerics::start(rr),
          allele("ref"), allele("alt"), sep = ":")
}

## The sample names by which combineRows() matches columns
.colKey <- function(obj) {
    key <- colnames(obj)
    if (is.null(key)) key <- as.character(colData(obj)$sample_id)
    key
}

## Put every object's rows or columns in the first object's order, by the
## keys that key_of() gives; an object whose keys are not the first
## object's is an error that names the argument and the keys that differ
.alignTo <- function(args, key_of, what, fun, hint = "") {
    ref <- key_of(args[[1L]])
    for (k in seq_along(args)[-1L]) {
        key <- key_of(args[[k]])
        if (identical(key, ref)) next
        if (anyDuplicated(ref) || anyDuplicated(key) ||
            !setequal(key, ref)) {
            lacks <- setdiff(ref, key)
            has <- setdiff(key, ref)
            show <- function(v) paste(utils::head(v, 5L), collapse = ", ")
            parts <- c(
                if (length(lacks))
                    sprintf("lacks %d of the first object's %s (%s)",
                            length(lacks), what, show(lacks)),
                if (length(has))
                    sprintf("has %d %s the first object lacks (%s)",
                            length(has), what, show(has)))
            if (!length(parts))
                parts <- paste("repeats some of its", what)
            stop(fun, " needs the same ", what, " in every object, but ",
                 "argument ", k, " ", paste(parts, collapse = " and "),
                 ".", hint, call. = FALSE)
        }
        idx <- match(ref, key)
        args[[k]] <- if (what == "samples") args[[k]][, idx]
                     else args[[k]][idx, ]
    }
    args
}

#' @keywords internal
#' @noRd
.mergeQcLogs <- function(whe_list) {
    logs_list <- lapply(whe_list, qcLog)
    ## Convert to data.frame for safe rbind, then back
    dfs <- lapply(logs_list, as.data.frame)
    merged <- do.call(base::rbind, dfs)
    if (nrow(merged) == 0L) {
        return(DataFrame(
            step = character(), parameter = character(),
            value = character(), n_flagged = integer(),
            n_passed = integer(), timestamp = character()))
    }
    DataFrame(merged)
}

# ============================================================
# Coercion methods -- interoperability
# ============================================================

#' Convert between WithinHostExperiment, data.frame, VRanges and VCF
#'
#' \code{as(x, "data.frame")} gives one row per site and sample with an
#' observed frequency; \code{as(x, "VRanges")} gives the same rows as a
#' \code{\link[VariantAnnotation]{VRanges}}. In the other direction,
#' \code{as(vr, "WithinHostExperiment")} builds the object from the
#' \code{refDepth}, \code{altDepth} and \code{totalDepth} of a
#' \code{VRanges}, one sample per \code{sampleNames} level, and
#' \code{as(vcf, "WithinHostExperiment")} from a VCF read with
#' \code{VariantAnnotation::readVcf()}, choosing the field layout from
#' the header as \code{\link{readWithinHost}(caller = "auto")} does and
#' taking every sample the file holds. Coercion to and from a
#' \code{VRanges} keeps the whole \code{seqinfo}: every recorded sequence
#' and its length, including sequences without a call, since the genome
#' length that \code{\link{calcDiversity}} divides by sums over all of
#' them. The VCF coercion keeps the lengths of the file's \code{##contig}
#' lines.
#'
#' @name coerce-WHE
#' @rdname coerce-WHE
#' @aliases coerce,WithinHostExperiment,data.frame-method
#'   coerce,WithinHostExperiment,VRanges-method
#'   coerce,VRanges,WithinHostExperiment-method
#'   coerce,VCF,WithinHostExperiment-method
#'
#' @return A \code{data.frame}, a
#'   \code{\link[VariantAnnotation]{VRanges}} or a
#'   \code{\link{WithinHostExperiment}}.
#'
#' @seealso \code{\link{readWithinHost}} for VCF files.
#'
#' @examples
#' data(example_whe)
#' df <- as(example_whe, "data.frame")
#' head(df)
#'
#' vr <- as(example_whe, "VRanges")
#' vr
#' as(vr, "WithinHostExperiment")
#'
#' vcf <- VariantAnnotation::readVcf(system.file("extdata",
#'     "test_donor.vcf", package = "WithinHostExperiment"))
#' as(vcf, "WithinHostExperiment")
NULL

setAs("WithinHostExperiment", "data.frame", function(from) {
    rr       <- rowRanges(from)
    freq_mat <- assay(from, "altFreq")
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

    n_sites   <- nrow(from)
    n_samples <- ncol(from)

    ref_vals <- if ("ref" %in% colnames(mcols(rr))) {
        mcols(rr)$ref
    } else {
        rep(NA_character_, n_sites)
    }
    alt_vals <- if ("alt" %in% colnames(mcols(rr))) {
        mcols(rr)$alt
    } else {
        rep(NA_character_, n_sites)
    }

    result <- data.frame(
        chrom     = rep(as.character(seqnames(rr)), times = n_samples),
        position  = rep(start(rr), times = n_samples),
        ref       = rep(ref_vals, times = n_samples),
        alt       = rep(alt_vals, times = n_samples),
        sample_id = rep(sample_ids, each = n_sites),
        altFreq   = as.numeric(freq_mat),
        stringsAsFactors = FALSE
    )

    for (aname in setdiff(assayNames(from), c("altFreq", "qcPass"))) {
        result[[aname]] <- as.numeric(assay(from, aname))
    }

    if ("qcPass" %in% assayNames(from)) {
        result$qcPass <- as.logical(assay(from, "qcPass"))
    }

    result <- result[!is.na(result$altFreq), , drop = FALSE]
    rownames(result) <- NULL
    result
})

setAs("WithinHostExperiment", "VRanges", function(from) {
    rr       <- rowRanges(from)
    freq_mat <- assay(from, "altFreq")
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

    ref_col <- if ("ref" %in% colnames(mcols(rr))) {
        mcols(rr)$ref
    } else {
        rep(NA_character_, length(rr))
    }
    alt_col <- if ("alt" %in% colnames(mcols(rr))) {
        mcols(rr)$alt
    } else {
        rep(NA_character_, length(rr))
    }

    depth_mat <- if ("totalDepth" %in% assayNames(from)) {
        assay(from, "totalDepth")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }
    alt_mat <- if ("altCount" %in% assayNames(from)) {
        assay(from, "altCount")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }
    ref_mat <- if ("refCount" %in% assayNames(from)) {
        assay(from, "refCount")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }

    n_sites   <- nrow(from)
    n_samples <- ncol(from)

    all_seqnames    <- rep(as.character(seqnames(rr)), times = n_samples)
    all_start       <- rep(start(rr), times = n_samples)
    all_ref         <- rep(ref_col, times = n_samples)
    all_alt         <- rep(alt_col, times = n_samples)
    all_totalDepth  <- as.integer(depth_mat)
    all_altDepth    <- as.integer(alt_mat)
    all_refDepth    <- as.integer(ref_mat)
    all_sampleNames <- rep(sample_ids, each = n_sites)
    all_freq        <- as.numeric(freq_mat)

    keep <- !is.na(all_freq)

    vr <- VRanges(
        seqnames    = factor(all_seqnames[keep],
                             levels = GenomeInfoDb::seqlevels(rr)),
        ranges      = IRanges(all_start[keep], width = 1L),
        ref         = all_ref[keep],
        alt         = all_alt[keep],
        totalDepth  = all_totalDepth[keep],
        altDepth    = all_altDepth[keep],
        refDepth    = all_refDepth[keep],
        sampleNames = all_sampleNames[keep]
    )
    ## keep every recorded sequence and its length, including those
    ## without a call, since the genome length sums over all of them
    GenomeInfoDb::seqinfo(vr) <- GenomeInfoDb::seqinfo(rr)
    vr
})


# ============================================================
# VRanges -> WithinHostExperiment coercion
# ============================================================

setAs("VRanges", "WithinHostExperiment", function(from) {
    ## Extract components
    sample_names <- as.character(
        VariantAnnotation::sampleNames(from))
    unique_samples <- unique(sample_names)
    n_samples <- length(unique_samples)

    ## Build site key for grouping
    chrom_vec <- as.character(GenomeInfoDb::seqnames(from))
    start_vec <- GenomicRanges::start(from)
    ref_vec <- as.character(VariantAnnotation::ref(from))
    alt_vec <- as.character(VariantAnnotation::alt(from))
    site_key <- paste(chrom_vec, start_vec, ref_vec, alt_vec,
                      sep = ":")
    unique_sites <- unique(site_key)
    n_sites <- length(unique_sites)

    ## One range per site, taken from its first record rather than parsed
    ## back out of the key, so a sequence name may contain ':'. The whole
    ## seqinfo is kept: every recorded sequence and its length.
    first <- match(unique_sites, site_key)
    gr <- GRanges(
        seqnames = factor(chrom_vec[first],
                          levels = GenomeInfoDb::seqlevels(from)),
        ranges = IRanges(start_vec[first], width = 1L),
        seqinfo = GenomeInfoDb::seqinfo(from)
    )
    mcols(gr)$ref <- ref_vec[first]
    mcols(gr)$alt <- alt_vec[first]

    ## Initialize matrices
    alt_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                         ncol = n_samples,
                         dimnames = list(NULL, unique_samples))
    ref_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                         ncol = n_samples,
                         dimnames = list(NULL, unique_samples))
    total_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                           ncol = n_samples,
                           dimnames = list(NULL, unique_samples))

    ## Fill matrices
    site_idx <- match(site_key, unique_sites)
    sample_idx <- match(sample_names, unique_samples)
    ad <- as.integer(VariantAnnotation::altDepth(from))
    rd <- as.integer(VariantAnnotation::refDepth(from))
    td <- as.integer(VariantAnnotation::totalDepth(from))
    for (k in seq_along(from)) {
        i <- site_idx[k]
        j <- sample_idx[k]
        alt_dp_mat[i, j] <- ad[k]
        ref_dp_mat[i, j] <- rd[k]
        total_dp_mat[i, j] <- td[k]
    }

    ## Compute frequency
    freq_mat <- alt_dp_mat / total_dp_mat

    WithinHostExperiment(
        assays = list(
            altFreq = freq_mat,
            totalDepth = total_dp_mat,
            altCount = alt_dp_mat,
            refCount = ref_dp_mat
        ),
        rowRanges = gr,
        colData = DataFrame(sample_id = unique_samples)
    )
})

# A misspelled argument absorbed by ... must not silently change an analysis.
.stop_on_unused_dots <- function(function_name, ...) {
    n_extra <- ...length()
    if (!n_extra) return(invisible(NULL))
    extra_names <- ...names()
    if (is.null(extra_names)) extra_names <- rep("", n_extra)
    extra_names[is.na(extra_names) | !nzchar(extra_names)] <- "<unnamed>"
    stop("Unused argument(s) in ", function_name, "(): ",
         paste(extra_names, collapse = ", "), ".", call. = FALSE)
}
