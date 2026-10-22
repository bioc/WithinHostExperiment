# R/import-vcf.R
# Import iSNV data from VCF files and VCF objects

#' @include methods.R annotate-codon.R
#' @importFrom VariantAnnotation readVcf expand header info geno samples
#'   ref alt
#' @importClassesFrom VariantAnnotation VCF
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqlengths seqlengths<- seqlevels
#'   seqlevelsInUse
#' @importFrom IRanges IRanges start end width
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
NULL

# ============================================================
# readWithinHost -- main entry point
# ============================================================

#' Import iSNV data from VCF files
#'
#' Reads one or more VCF files produced by variant callers (iVar
#' through a VCF converter, LoFreq, Freebayes, or others) with
#' \code{VariantAnnotation::readVcf()} and constructs a
#' \code{\link{WithinHostExperiment}}. A VCF already read into R is
#' converted with \code{as(vcf, "WithinHostExperiment")}, and a
#' \code{VRanges} with \code{as(vr, "WithinHostExperiment")}; see
#' \code{\link{coerce-WHE}}.
#'
#' Each VCF file corresponds to one sample. The \code{colData}
#' argument provides sample metadata and must contain a
#' \code{sample_id} column whose length matches the number of files. A
#' file holding more than one sample is an error; read it with
#' \code{readVcf()} and convert the object instead.
#'
#' \code{caller = "ivar"} refers to a VCF carrying iVar's per-sample
#' fields (\code{REF_DP}, \code{ALT_DP}, \code{ALT_FREQ},
#' \code{ALT_QUAL}). \code{ivar variants} itself writes a TSV, not a
#' VCF, so such files come from a converter, for instance
#' \code{ivar_variants_to_vcf.py} in the nf-core/viralrecon pipeline.
#' An iVar TSV is read directly by
#' \code{\link{readWithinHostTable}(format = "ivar")}.
#'
#' @param vcfFiles Character vector. Paths to VCF files (one per sample).
#' @param colData \code{DataFrame} or \code{data.frame}. Sample
#'   metadata. Must contain \code{sample_id}. May also contain
#'   \code{host_id}, \code{role}, \code{pair_id},
#'   \code{replicate_group}, \code{replicate_id}, etc.
#' @param caller Character scalar. Which field layout the VCF uses:
#'   \code{"ivar"}, \code{"lofreq"}, \code{"freebayes"}, \code{"auto"}
#'   (default) or \code{"generic"}. See Details.
#' @param genome Optional \code{DNAStringSet} holding the reference
#'   sequence, as returned by \code{Biostrings::readDNAStringSet()}. The
#'   length of each sequence is recorded in
#'   \code{seqlengths(rowRanges())}. It must name the sequences the calls
#'   lie on and be long enough to hold them; see Details.
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @details
#' \strong{Caller.} With \code{caller = "auto"} the layout is chosen
#' from the VCF header: a header naming iVar, LoFreq or Freebayes (as
#' the \code{##source} line of each does) selects that layout, and any
#' other gives \code{"generic"}. The layouts read: iVar, FORMAT
#' \code{ALT_FREQ}, \code{ALT_DP} and \code{REF_DP}, with depth their
#' sum (INFO \code{DP} where they are missing); LoFreq, INFO \code{AF}
#' and \code{DP}, with counts derived from them; Freebayes, FORMAT
#' \code{AO}, \code{RO} and \code{DP} (INFO \code{DP} where FORMAT has
#' none), with frequency \code{AO / DP}; generic, the first of those
#' field names that the file carries.
#'
#' \strong{Several alternative alleles.} A record with several
#' alternative alleles is split into one site per allele with
#' \code{VariantAnnotation::expand()}. A per-allele value is split
#' correctly only when the header declares it \code{Number=A}, as the
#' VCF specification requires; \code{readVcf()} keeps just the first
#' value of a field declared otherwise. When such a record carries one
#' of the per-allele fields read here (INFO \code{AF}; FORMAT
#' \code{AO}, \code{AF}, \code{ALT_DP}, \code{ALT_FREQ},
#' \code{ALT_QUAL}) under another \code{Number}, the function stops
#' rather than give the second allele the first allele's value. Correct
#' the header, or split the records first (for example with
#' \code{bcftools norm -m-}).
#'
#' \strong{Sequence lengths.} Lengths are recorded in
#' \code{seqlengths(rowRanges())} from the VCF header's \code{##contig}
#' lines and from \code{genome}; two sources that disagree are an
#' error. Functions that need the genome length read it from there, so
#' \code{genome} is checked before anything is recorded: it is an error
#' when it names none of the sequences the calls lie on, or when a call
#' lies past the end of its sequence, since either way it is not the
#' reference the calls were made against. A sequence the calls lie on
#' that \code{genome} does not name is left without a length, with a
#' warning.
#'
#' When technical replicates are not available, Roder et al. (2023,
#' \emph{mBio} 14:e01046-23) recommend using a combination of
#' multiple variant callers with stringent cutoffs to reduce false
#' positives. To implement multi-caller consensus, import VCFs from
#' each caller separately and intersect the resulting objects using
#' standard \code{GRanges} overlap operations.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' whe
#'
#' ## The same file read into R first, then converted
#' v <- VariantAnnotation::readVcf(vcf)
#' as(v, "WithinHostExperiment")
readWithinHost <- function(vcfFiles, colData, caller = "auto",
                           genome = NULL) {
    ## ---- Validate inputs ----
    if (!is.character(vcfFiles) || length(vcfFiles) < 1L)
        stop("'vcfFiles' must be a non-empty character vector.")
    missing_f <- vcfFiles[!file.exists(vcfFiles)]
    if (length(missing_f) > 0L) {
        stop("VCF file(s) not found: ", paste(missing_f, collapse = ", "))
    }

    if (is.data.frame(colData) && !is(colData, "DataFrame")) {
        colData <- DataFrame(colData)
    }
    if (!is(colData, "DataFrame"))
        stop("'colData' must be a DataFrame (S4Vectors::DataFrame).")

    if (!"sample_id" %in% colnames(colData)) {
        stop("'colData' must contain a 'sample_id' column. ",
             "Found columns: ",
             paste(colnames(colData), collapse = ", "))
    }
    if (nrow(colData) != length(vcfFiles)) {
        stop("Length of 'vcfFiles' (", length(vcfFiles),
             ") must match nrow(colData) (", nrow(colData), ")")
    }
    if (is.character(genome))
        stop("'genome' must be a DNAStringSet, not a file path. Read the ",
             "FASTA file first, for example genome <- ",
             "Biostrings::readDNAStringSet(\"", genome[1L], "\").",
             call. = FALSE)
    if (!is.null(genome) && !is(genome, "DNAStringSet"))
        stop("'genome' must be a DNAStringSet or NULL.")

    valid_callers <- c("ivar", "lofreq", "freebayes", "auto", "generic")
    caller <- match.arg(caller, valid_callers)

    ## ---- Read each VCF ----
    vcfs <- lapply(vcfFiles, function(f) readVcf(f))
    n_samp <- vapply(vcfs, function(v) length(samples(header(v))),
                     integer(1))
    if (any(n_samp > 1L))
        stop("readWithinHost() reads one sample per file, but ",
             paste0(basename(vcfFiles[n_samp > 1L]), " holds ",
                    n_samp[n_samp > 1L], " samples", collapse = "; "),
             ". Read a multi-sample file with VariantAnnotation::readVcf() ",
             "and convert it with as(vcf, \"WithinHostExperiment\").",
             call. = FALSE)

    if (caller == "auto") {
        caller <- .detectCaller(header(vcfs[[1L]]))
        message("Auto-detected caller: ", caller)
    }

    parsed <- lapply(seq_along(vcfs), function(i) {
        .checkPerAlleleNumbers(vcfs[[i]], basename(vcfFiles[i]))
        .parseVcfObject(vcfs[[i]], caller)[[1L]]
    })
    lens <- .combineSeqlengths(lapply(vcfs, seqlengths))
    whe <- .mergeToWHE(parsed, colData, callerName = caller,
                       seqlengths = lens)
    if (!is.null(genome)) {
        genome <- .prepareReference(genome)
        .checkGenome(rowRanges(whe), genome)
        SummarizedExperiment::rowRanges(whe) <-
            .fillSeqlengths(rowRanges(whe), genome)
    }
    whe
}

## Stop when 'genome' names none of the sequences the calls lie on, or
## when a call lies past the end of its sequence: either way it is not
## the reference the calls were made against, and recording its lengths
## would give every per-site statistic a wrong denominator. A sequence
## it does not name is left without a length, with a warning.
.checkGenome <- function(rr, genome) {
    v <- seqlevelsInUse(rr)
    if (length(v) == 0L) return(invisible(TRUE))
    g <- names(genome)
    if (length(intersect(v, g)) == 0L)
        stop("'genome' names none of the sequences the calls lie on, so ",
             "no length could be recorded.\n",
             "  variant seqlevels: ", paste(v, collapse = ", "), "\n",
             "  'genome' sequences: ", paste(g, collapse = ", "), "\n",
             "Rename one of them to match, for example with ",
             "GenomeInfoDb::renameSeqlevels().", call. = FALSE)
    unnamed <- setdiff(v, g)
    if (length(unnamed))
        warning("'genome' has no sequence named ",
                paste(unnamed, collapse = ", "),
                "; no length is recorded for it.", call. = FALSE)
    chrom <- as.character(seqnames(rr))
    len <- width(genome)[match(chrom, g)]
    beyond <- !is.na(len) & end(rr) > len
    if (any(beyond)) {
        k <- which(beyond)[1L]
        stop(sum(beyond), " call(s) lie past the end of their 'genome' ",
             "sequence, the first at ", chrom[k], ":", start(rr)[k],
             " on a sequence of ", len[k], " bases, so 'genome' is not ",
             "the reference the calls were made against.", call. = FALSE)
    }
    invisible(TRUE)
}


# ============================================================
# Internal: caller detection
# ============================================================

## The header's metadata lines, names and values, as one string: the
## ##source line normally names the caller. Takes a VCFHeader or a path.
#' @keywords internal
.detectCaller <- function(hdr) {
    if (is.character(hdr)) hdr <- VariantAnnotation::scanVcfHeader(hdr)
    m <- VariantAnnotation::meta(hdr)
    fields <- unlist(lapply(m, function(d)
        c(rownames(d), unlist(lapply(as.list(d), as.character)))),
        use.names = FALSE)
    text <- paste(c(names(m), fields), collapse = "\n")
    if (grepl("iVar", text, ignore.case = TRUE)) return("ivar")
    if (grepl("lofreq", text, ignore.case = TRUE)) return("lofreq")
    if (grepl("freebayes", text, ignore.case = TRUE)) return("freebayes")
    "generic"
}


# ============================================================
# Internal: reading values from a VCF object
# ============================================================

.PER_ALLELE_INFO <- "AF"
.PER_ALLELE_FORMAT <- c("AO", "AF", "ALT_DP", "ALT_FREQ", "ALT_QUAL")

## readVcf() keeps only the first value of a field whose header does not
## declare it per allele, so a multi-allelic record would give every allele
## the first allele's value. Stop instead.
#' @keywords internal
.checkPerAlleleNumbers <- function(vcf, label = "the VCF") {
    multi <- lengths(alt(vcf)) > 1L
    if (!any(multi)) return(invisible(TRUE))
    hdr <- header(vcf)
    num <- function(d, ids) {
        ids <- intersect(ids, rownames(d))
        stats::setNames(as.character(d[ids, "Number"]), ids)
    }
    declared <- c(num(info(hdr), .PER_ALLELE_INFO),
                  num(geno(hdr), .PER_ALLELE_FORMAT))
    bad <- declared[!declared %in% c("A", "R")]
    if (length(bad))
        stop(label, " has ", sum(multi), " record(s) with several ",
             "alternative alleles, but its header declares ",
             paste0(names(bad), " as Number=", bad, collapse = ", "),
             " rather than Number=A, so the value of each allele cannot ",
             "be read. Correct the header, or split the records first, ",
             "for example with 'bcftools norm -m-'.", call. = FALSE)
    invisible(TRUE)
}

## One value per site from an INFO column, or NA when absent.
#' @keywords internal
.infoValue <- function(e, tag, n) {
    d <- info(e)
    if (!tag %in% colnames(d)) return(rep(NA_real_, n))
    v <- d[[tag]]
    if (is(v, "List") || is.list(v))
        v <- vapply(as.list(v), function(z)
            if (length(z)) as.numeric(z[[1L]]) else NA_real_, numeric(1))
    as.numeric(v)
}

## One value per site for sample j from a FORMAT field, or NA when absent.
#' @keywords internal
.genoValue <- function(e, field, j, n) {
    g <- geno(e)
    if (!field %in% names(g)) return(rep(NA_real_, n))
    m <- g[[field]]
    v <- if (length(dim(m)) == 3L) m[, j, 1L] else m[, j]
    if (is.list(v))
        v <- vapply(v, function(z)
            if (length(z)) as.character(z[[1L]]) else NA_character_,
            character(1))
    if (is.numeric(v)) return(as.numeric(v))
    ## a missing value is "." in a VCF; anything else that is not a
    ## number is NA as well
    v <- as.character(v)
    number <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
    numeric_text <- grepl(number, v)
    out <- rep(NA_real_, length(v))
    out[numeric_text] <- as.numeric(v[numeric_text])
    out
}

## Per-sample site tables in the layout .mergeToWHE() takes: one element
## per sample, or one element for a VCF without FORMAT columns.
#' @keywords internal
.parseVcfObject <- function(vcf, caller) {
    e <- expand(vcf)
    n <- nrow(e)
    rr <- SummarizedExperiment::rowRanges(e)
    sites <- list(
        chrom = as.character(seqnames(rr)),
        pos   = as.integer(start(rr)),
        ref   = as.character(ref(e)),
        alt   = as.character(alt(e)))
    as_int <- function(v) as.integer(round(v))
    info_dp <- as_int(.infoValue(e, "DP", n))
    info_only <- function() {
        alt_freq <- .infoValue(e, "AF", n)
        total_dp <- info_dp
        alt_dp <- as_int(total_dp * alt_freq)
        list(altFreq = alt_freq, totalDepth = total_dp,
             refCount = total_dp - alt_dp, altCount = alt_dp)
    }
    n_samples <- length(samples(header(vcf)))
    if (caller == "lofreq" || n_samples == 0L)
        return(list(c(sites, info_only())))

    lapply(seq_len(n_samples), function(j) {
        fmt <- function(field) .genoValue(e, field, j, n)
        first_of <- function(fields) {
            for (f in fields) {
                v <- fmt(f)
                if (!all(is.na(v))) return(v)
            }
            rep(NA_real_, n)
        }
        vals <- if (caller == "ivar") {
            alt_dp <- as_int(fmt("ALT_DP"))
            ref_dp <- as_int(fmt("REF_DP"))
            total <- alt_dp + ref_dp
            list(altFreq = fmt("ALT_FREQ"),
                 totalDepth = ifelse(is.na(total), info_dp, total),
                 refCount = ref_dp, altCount = alt_dp)
        } else if (caller == "freebayes") {
            alt_dp <- as_int(fmt("AO"))
            ref_dp <- as_int(fmt("RO"))
            total <- as_int(fmt("DP"))
            if (all(is.na(total))) total <- info_dp
            list(altFreq = ifelse(!is.na(alt_dp) & !is.na(total) & total > 0L,
                                  alt_dp / total, NA_real_),
                 totalDepth = total, refCount = ref_dp, altCount = alt_dp)
        } else {
            alt_dp <- as_int(first_of(c("ALT_DP", "AO")))
            ref_dp <- as_int(first_of(c("REF_DP", "RO")))
            total <- as_int(fmt("DP"))
            if (all(is.na(total)))
                total <- ifelse(!is.na(alt_dp) & !is.na(ref_dp),
                                alt_dp + ref_dp, info_dp)
            list(altFreq = first_of(c("ALT_FREQ", "AF")),
                 totalDepth = total, refCount = ref_dp, altCount = alt_dp)
        }
        c(sites, vals)
    })
}

## Merge named sequence lengths from several sources; NA means unknown,
## and two known lengths that differ are an error.
#' @keywords internal
.combineSeqlengths <- function(lens_list) {
    out <- integer()
    for (l in lens_list) {
        for (nm in names(l)) {
            v <- l[[nm]]
            if (is.na(v)) {
                if (!nm %in% names(out)) out[nm] <- NA_integer_
            } else if (!nm %in% names(out) || is.na(out[[nm]])) {
                out[nm] <- as.integer(v)
            } else if (out[[nm]] != v) {
                stop("The VCF headers give different lengths for ", nm,
                     ": ", out[[nm]], " and ", v, ".", call. = FALSE)
            }
        }
    }
    out
}


# ============================================================
# VCF object -> WithinHostExperiment coercion
# ============================================================

setAs("VCF", "WithinHostExperiment", function(from) {
    caller <- .detectCaller(header(from))
    .checkPerAlleleNumbers(from)
    parsed <- .parseVcfObject(from, caller)
    ids <- samples(header(from))
    if (!length(ids)) ids <- "sample1"
    .mergeToWHE(parsed, DataFrame(sample_id = ids), callerName = caller,
                seqlengths = .combineSeqlengths(list(seqlengths(from))))
})


# ============================================================
# Internal: merge multiple samples -> WithinHostExperiment
# ============================================================

#' @keywords internal
.mergeToWHE <- function(parsed_list, colData, callerName, metadata = list(),
                        seqlengths = NULL) {
    sample_ids <- as.character(colData$sample_id)
    n_samples  <- length(parsed_list)

    ## 1. Collect all unique positions across samples
    all_positions <- unique(do.call(rbind, lapply(parsed_list, function(p) {
        if (length(p$chrom) == 0L) {
            return(data.frame(chrom = character(), pos = integer(),
                              ref = character(), alt = character(),
                              stringsAsFactors = FALSE))
        }
        data.frame(chrom = p$chrom, pos = p$pos, ref = p$ref, alt = p$alt,
                   stringsAsFactors = FALSE)
    })))

    if (nrow(all_positions) == 0L) {
        ## No variants found - return empty WHE
        gr <- GRanges()
        empty_num <- matrix(nrow = 0, ncol = n_samples,
                            dimnames = list(NULL, sample_ids))
        empty_int <- matrix(integer(0), nrow = 0, ncol = n_samples,
                            dimnames = list(NULL, sample_ids))
        return(.makeWHE(
            rowRanges = gr,
            assays = list(altFreq = empty_num, totalDepth = empty_int,
                          refCount = empty_int, altCount = empty_int),
            colData = colData, callerName = callerName, metadata = metadata
        ))
    }

    ## Sort by chrom + pos
    all_positions <- all_positions[order(all_positions$chrom,
                                        all_positions$pos), ]
    ## Remove exact duplicates (same chrom+pos+ref+alt)
    dup_key <- paste(all_positions$chrom, all_positions$pos,
                     all_positions$ref, all_positions$alt, sep = ":")
    all_positions <- all_positions[!duplicated(dup_key), ]
    n_sites <- nrow(all_positions)

    ## 2. Build GRanges
    gr <- GRanges(
        seqnames = all_positions$chrom,
        ranges   = IRanges(all_positions$pos, width = 1L)
    )
    S4Vectors::mcols(gr)$ref <- all_positions$ref
    S4Vectors::mcols(gr)$alt <- all_positions$alt
    if (length(seqlengths)) {
        known <- seqlengths[intersect(names(seqlengths), seqlevels(gr))]
        if (length(known)) {
            lens <- seqlengths(gr)
            lens[names(known)] <- known
            seqlengths(gr) <- lens
        }
    }

    ## 3. Build assay matrices
    freq_mat  <- matrix(NA_real_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    depth_mat <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    ref_mat   <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    alt_mat   <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))

    master_key <- paste(all_positions$chrom, all_positions$pos,
                        all_positions$ref, all_positions$alt, sep = ":")

    for (i in seq_len(n_samples)) {
        p <- parsed_list[[i]]
        if (length(p$chrom) == 0L) next
        sample_key <- paste(p$chrom, p$pos, p$ref, p$alt, sep = ":")
        idx   <- match(sample_key, master_key)
        valid <- !is.na(idx)
        freq_mat[idx[valid], i]  <- p$altFreq[valid]
        depth_mat[idx[valid], i] <- p$totalDepth[valid]
        ref_mat[idx[valid], i]   <- p$refCount[valid]
        alt_mat[idx[valid], i]   <- p$altCount[valid]
    }

    ## 4. Construct WHE
    .makeWHE(
        rowRanges = gr,
        assays = list(altFreq    = freq_mat,
                      totalDepth = depth_mat,
                      refCount   = ref_mat,
                      altCount   = alt_mat),
        colData    = colData,
        callerName = callerName,
        metadata   = metadata
    )
}
