# R/import-table.R
# Import iSNV data from tabular files (iVar TSV, generic CSV)

#' @include methods.R
#' @importFrom utils read.csv read.delim
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
NULL

#' Import iSNV data from tabular files
#'
#' Reads pre-computed iSNV data from TSV or CSV files and constructs
#' a \code{\link{WithinHostExperiment}}.
#'
#' @param files Character vector. Paths to files (one per sample).
#' @param format Character. \code{"ivar"} for iVar TSV output (columns:
#'   REGION, POS, REF, ALT, REF_DP, ALT_DP, ALT_FREQ, TOTAL_DP, ...),
#'   or \code{"generic"} for CSV with columns: chrom, pos, ref, alt,
#'   aaf, depth.
#' @param colData \code{DataFrame} (optional). Sample metadata with
#'   \code{sample_id}. If \code{NULL}, sample IDs are derived from
#'   file basenames.
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @details
#' Only the columns named under \code{format} are read.
#'
#' iVar's \code{PASS} column, the outcome of its own test of each call,
#' is not read, so a row with \code{PASS} \code{FALSE} is imported as a
#' call, and \code{\link{flagISNV}} judges it by its depth, frequency and
#' read-count thresholds alone. To drop such rows, remove them from the
#' file first, reading every column as text so that the file is written
#' back unchanged:
#' \preformatted{
#' calls <- read.delim(path, colClasses = "character")
#' passed <- tempfile(fileext = ".tsv")
#' write.table(calls[calls$PASS == "TRUE", ], passed, sep = "\\t",
#'     quote = FALSE, row.names = FALSE)
#' whe <- readWithinHostTable(passed, format = "ivar")
#' }
#'
#' iVar's per-strand counts (\code{REF_RV}, \code{ALT_RV}) are not
#' imported. The strand-bias step of \code{flagISNV()} runs only when the
#' assays \code{fwdAltCount} and \code{revAltCount} are present, so on an
#' imported object it is skipped and leaves no \code{qcLog} row, although
#' \code{qcFilters()} records the filter's \code{maxStrandBias};
#' \code{\link{flagStrandBiasSCI}} stops for want of those assays.
#'
#' Insertions and deletions are imported as rows of width 1 at their
#' position, with the caller's allele notation (iVar writes \code{+TTG}
#' for an insertion and \code{-T} for a deletion), and nothing tells
#' them apart from substitutions: they pass \code{flagISNV()} like any
#' other call and enter the diversity statistics. A VCF read by
#' \code{\link{readWithinHost}} writes the same one-base deletion as, for
#' instance, \code{AT} to \code{A}, so the two routes do not match it by
#' allele. To keep substitutions only, subset to rows whose \code{ref}
#' and \code{alt} are one base each:
#' \preformatted{
#' rr <- SummarizedExperiment::rowRanges(whe)
#' whe <- whe[nchar(rr$ref) == 1L & nchar(rr$alt) == 1L, ]
#' }
#'
#' @export
#' @examples
#' tsv <- system.file("extdata", "test_donor_ivar.tsv",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHostTable(tsv, format = "ivar")
#' whe
#'
#' ## Only the rows iVar passed
#' calls <- read.delim(tsv, colClasses = "character")
#' passed <- tempfile(fileext = ".tsv")
#' write.table(calls[calls$PASS == "TRUE", ], passed, sep = "\t",
#'     quote = FALSE, row.names = FALSE)
#' readWithinHostTable(passed, format = "ivar")
readWithinHostTable <- function(files, format = c("ivar", "generic"),
                                colData = NULL) {
    format <- match.arg(format)
    if (!is.character(files) || length(files) < 1L)
        stop("'files' must be a non-empty character vector.")
    missing_f <- files[!file.exists(files)]
    if (length(missing_f) > 0L) {
        stop("File(s) not found: ", paste(missing_f, collapse = ", "))
    }

    ## Auto-generate colData if not provided
    if (is.null(colData)) {
        ids <- tools::file_path_sans_ext(basename(files))
        colData <- DataFrame(sample_id = ids)
    } else {
        if (is.data.frame(colData) && !is(colData, "DataFrame")) {
            colData <- DataFrame(colData)
        }
        if (!"sample_id" %in% colnames(colData)) {
            stop("'colData' must contain a 'sample_id' column. ",
                 "Found columns: ",
                 paste(colnames(colData), collapse = ", "))
        }
        if (nrow(colData) != length(files)) {
            stop("Length of 'files' (", length(files),
                 ") must match nrow(colData) (", nrow(colData), ")")
        }
    }

    ## Parse each file
    parsed <- lapply(files, function(f) {
        if (format == "ivar") {
            .parseIvarTsv(f)
        } else {
            .parseGenericCsv(f)
        }
    })

    .mergeToWHE(parsed, colData, callerName = format, metadata = list())
}


#' @keywords internal
.parseIvarTsv <- function(file) {
    ## read.delim() would turn a column that holds only "T" into TRUE
    df <- read.delim(file, stringsAsFactors = FALSE,
                     colClasses = c(REGION = "character", REF = "character",
                                    ALT = "character"))
    ## Expected columns: REGION, POS, REF, ALT, REF_DP, ALT_DP,
    ## ALT_FREQ, TOTAL_DP
    required <- c("REGION", "POS", "REF", "ALT")
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0L) {
        stop("iVar TSV missing columns: ",
             paste(missing_cols, collapse = ", "))
    }

    list(
        chrom      = df$REGION,
        pos        = as.integer(df$POS),
        ref        = df$REF,
        alt        = df$ALT,
        altFreq    = if ("ALT_FREQ" %in% colnames(df))
                         as.numeric(df$ALT_FREQ) else NA_real_,
        totalDepth = if ("TOTAL_DP" %in% colnames(df))
                         as.integer(df$TOTAL_DP) else NA_integer_,
        refCount   = if ("REF_DP" %in% colnames(df))
                         as.integer(df$REF_DP) else NA_integer_,
        altCount   = if ("ALT_DP" %in% colnames(df))
                         as.integer(df$ALT_DP) else NA_integer_
    )
}


#' @keywords internal
.parseGenericCsv <- function(file) {
    ## read.csv() would turn a column that holds only "T" into TRUE
    df <- read.csv(file, stringsAsFactors = FALSE,
                   colClasses = c(chrom = "character", ref = "character",
                                  alt = "character"))
    required <- c("chrom", "pos", "ref", "alt", "aaf", "depth")
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0L) {
        stop("Generic CSV missing columns: ",
             paste(missing_cols, collapse = ", "))
    }

    alt_dp <- as.integer(round(df$depth * df$aaf))
    ref_dp <- as.integer(df$depth) - alt_dp

    list(
        chrom      = df$chrom,
        pos        = as.integer(df$pos),
        ref        = df$ref,
        alt        = df$alt,
        altFreq    = as.numeric(df$aaf),
        totalDepth = as.integer(df$depth),
        refCount   = ref_dp,
        altCount   = alt_dp
    )
}
