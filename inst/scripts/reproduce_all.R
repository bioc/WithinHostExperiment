#!/usr/bin/env Rscript
## reproduce_all.R
##
## Rebuilds every derived file in the repository, from the data in
## inst/scripts/real_data to the README figures, and checks each one
## against the copy that is committed. Nothing is typed by hand at any
## stage: each derived file is written by one script listed below.
##
## Run from the package root of a clone:
##   Rscript inst/scripts/reproduce_all.R             stages 1-10
##   Rscript inst/scripts/reproduce_all.R --download  also stage 0
##
## Stages, each run as its own R process with its log kept:
##   0. prepare_farjo_metadata.R   (--download only; needs network and
##                                 readxl) samples_432870.csv from the
##                                 Farjo et al. repository at a pinned
##                                 commit
##   1. inputs                     every input file against the checksums
##                                 committed in inst/scripts/readme_values/
##   2. install                    this source into a temporary library
##   3. prepare_lauring_data.R     Bendall et al. excerpts: real_data/MHM*
##                                 and the four HH6 files in inst/extdata
##   4. create_test_data.R         the twelve synthetic fixtures (seed 42),
##                                 then the three fixtures without a
##                                 generator against the checksums that
##                                 make_extdata.R records
##   5. save_example_data.R        data/example_whe.rda
##   6. tests                      tests/testthat against the installed
##                                 package
##   7. generate_readme_figures.R  man/figures/*.png, the PDF copies and
##                                 every number the captions quote
##   8. fill_readme.R              README.md from its template
##   9. verify_claims.R            the functional checks, against the
##                                 installed package
##
## A stage passes when its script exits without error and every file it
## writes is byte-identical to what was there before. data/example_whe.rda
## is compared by its content (assays, rowRanges, colData, metadata, qcLog,
## qcFilters), since the file also records the R version and the packages
## that defined its classes. sessionInfo.txt records the run's own packages
## and is reported, not judged. The vector PDFs are not tracked and are
## only listed. The run exits with status 1 if any stage fails and leaves
## what it wrote in place for inspection; git restores the committed copies.

args <- commandArgs(trailingOnly = TRUE)
download <- "--download" %in% args
unknown <- setdiff(args, "--download")
if (length(unknown))
    stop("Unknown argument(s): ", paste(unknown, collapse = " "))

if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] !=
        "WithinHostExperiment")
    stop("Run this script from the WithinHostExperiment package root.")

rscript <- file.path(R.home("bin"), "Rscript")
rcmd <- file.path(R.home("bin"), "R")
## Logs go beside the session's temporary directory, which R removes on exit
log_dir <- file.path(dirname(tempdir()), "WithinHostExperiment_reproduce")
unlink(log_dir, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
lib_dir <- file.path(log_dir, "library")
dir.create(lib_dir, showWarnings = FALSE)
lib_env <- c(R_LIBS = paste(c(lib_dir, .libPaths()),
                            collapse = .Platform$path.sep))

DATA <- file.path("inst", "scripts", "real_data")
EXT <- file.path("inst", "extdata")
OUT <- file.path("inst", "scripts", "readme_values")

## The checksum of a text file is taken with carriage returns removed, so
## that it is the checksum of the committed file on any platform: a Windows
## checkout with core.autocrlf adds a carriage return to every line.
md5_text <- function(files) {
    vapply(files, function(f) {
        if (!file.exists(f)) return(NA_character_)
        b <- readBin(f, "raw", file.size(f))
        tf <- tempfile()
        on.exit(unlink(tf))
        writeBin(b[b != as.raw(13L)], tf)
        unname(tools::md5sum(tf))
    }, character(1), USE.NAMES = FALSE)
}

## Byte checksums, for comparing a file with itself before and after a
## stage in the same checkout
md5 <- function(files) {
    v <- rep(NA_character_, length(files))
    ok <- file.exists(files)
    v[ok] <- unname(tools::md5sum(files[ok]))
    stats::setNames(v, files)
}

results <- data.frame(stage = character(), status = character(),
                      detail = character(), stringsAsFactors = FALSE)
note <- function(stage, status, detail = "") {
    results[nrow(results) + 1L, ] <<- list(stage, status, detail)
    cat(sprintf("[%s] %-27s %s\n", status, stage, detail))
}

## Run an R script in a child process. system2(env =) does not reach the
## child on Windows, so variables are set in this process around the call.
run_r <- function(script, log, env = character()) {
    old <- Sys.getenv(names(env), unset = NA)
    if (length(env)) do.call(Sys.setenv, as.list(env))
    on.exit({
        for (nm in names(old)) {
            if (is.na(old[[nm]])) Sys.unsetenv(nm)
            else do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
        }
    })
    system2(rscript, shQuote(script), stdout = log, stderr = log)
}

## Run one stage; compare the files it writes before and after. With
## announce = FALSE the caller reports the outcome itself.
run_stage <- function(stage, script, watch = character(), env = character(),
                      report_only = character(), announce = TRUE) {
    before <- md5(c(watch, report_only))
    log <- file.path(log_dir, paste0(gsub("[^A-Za-z0-9]+", "_", stage),
                                     ".log"))
    status <- run_r(script, log, env)
    after <- md5(c(watch, report_only))
    if (!identical(as.integer(status), 0L)) {
        note(stage, "FAIL", paste0("exit ", status, "; log ", log))
        return(invisible(FALSE))
    }
    differs <- function(f) is.na(before[f]) | is.na(after[f]) |
        before[f] != after[f]
    changed <- watch[differs(watch)]
    moved <- report_only[differs(report_only)]
    if (length(changed)) {
        note(stage, "FAIL", paste("changed:", paste(changed, collapse = ", ")))
        return(invisible(FALSE))
    }
    if (announce)
        note(stage, "PASS", paste0(
            if (length(watch)) paste(length(watch), "file(s) unchanged")
            else "completed",
            if (length(moved)) paste0("; environment record updated: ",
                                      paste(moved, collapse = ", "))))
    invisible(TRUE)
}

## A small script written to the log directory and run like the others.
helper <- function(name, lines) {
    f <- file.path(log_dir, name)
    writeLines(lines, f)
    f
}

cat("WithinHostExperiment: rebuild from raw data\n")
cat(R.version.string, "| Bioconductor",
    tryCatch(as.character(BiocManager::version()),
             error = function(e) "unknown"),
    "\nlogs:", log_dir, "\n\n")

## ---- 0. Farjo metadata from the upstream repository -------------------
farjo_meta <- file.path(DATA, "farjo_longitudinal", "samples_432870.csv")
if (!download) {
    note("0 prepare_farjo_metadata", "SKIP", "pass --download to run it")
} else if (!requireNamespace("readxl", quietly = TRUE)) {
    note("0 prepare_farjo_metadata", "SKIP", "readxl is not installed")
} else {
    run_stage("0 prepare_farjo_metadata",
              file.path("inst", "scripts", "prepare_farjo_metadata.R"),
              watch = farjo_meta)
}

## ---- 1. Inputs against the committed checksums --------------------------
recorded <- utils::read.csv(file.path(OUT, "input_md5.csv"),
                            stringsAsFactors = FALSE)
now <- md5_text(recorded$file)
bad <- recorded$file[is.na(now) | now != recorded$md5]
if (length(bad)) {
    note("1 inputs", "FAIL", paste("differ from input_md5.csv:",
                                   paste(bad, collapse = ", ")))
} else {
    note("1 inputs", "PASS", paste(nrow(recorded),
                                   "input file(s) match input_md5.csv"))
}

## ---- 2. Install this source -------------------------------------------
install_log <- file.path(log_dir, "2_install.log")
st <- system2(rcmd, c("CMD", "INSTALL", "--no-test-load", "-l",
                      shQuote(lib_dir), "."),
              stdout = install_log, stderr = install_log)
installed <- identical(as.integer(st), 0L)
note("2 install", if (installed) "PASS" else "FAIL",
     if (installed) lib_dir else paste0("exit ", st, "; log ", install_log))

## ---- 3. Bendall et al. excerpts ---------------------------------------
lauring_out <- c(
    list.files(DATA, pattern = "^MHM.*\\.tsv$", full.names = TRUE),
    file.path(EXT, paste0(rep(c("MHM1143", "MHM1184"), each = 2),
                          c("_rep1.tsv", "_rep2.tsv"))))
run_stage("3 prepare_lauring_data",
          file.path("inst", "scripts", "prepare_lauring_data.R"),
          watch = lauring_out)

## ---- 4. Synthetic fixtures --------------------------------------------
fixtures <- file.path(EXT, c("test_ref.fa", "test_donor.vcf",
    "test_recipient.vcf", "test_independent.vcf", "test_donor_rep1.vcf",
    "test_donor_rep2.vcf", "test_donor_ivar.tsv", "test_donor_lofreq.vcf",
    "test_donor_freebayes.vcf", "test_generic.csv", "test_coldata.csv",
    "test_primers.bed"))
run_stage("4 create_test_data",
          file.path("inst", "scripts", "create_test_data.R"),
          watch = fixtures)

## The three fixtures without a generator. Their checksums are recorded
## once, in make_extdata.R, and read from it here.
extdata_exprs <- parse(file.path("inst", "scripts", "make_extdata.R"))
is_unchanged <- vapply(extdata_exprs, function(e)
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("unchanged")), logical(1))
unchanged <- eval(extdata_exprs[[which(is_unchanged)]][[3L]])
now <- md5_text(file.path(EXT, names(unchanged)))
moved <- names(unchanged)[is.na(now) | now != unname(unchanged)]
note("4 fixtures, no generator", if (length(moved)) "FAIL" else "PASS",
     if (length(moved)) paste("changed:", paste(moved, collapse = ", "))
     else paste(length(unchanged), "file(s) match make_extdata.R"))

## ---- 5. Example dataset, compared by content --------------------------
rda <- file.path("data", "example_whe.rda")
rda_before <- file.path(log_dir, "example_whe_before.rda")
invisible(file.copy(rda, rda_before, overwrite = TRUE))
if (run_stage("5 save_example_data",
              file.path("inst", "scripts", "save_example_data.R"),
              announce = FALSE)) {
    cmp <- helper("compare_example_whe.R", c(
        "suppressPackageStartupMessages({",
        "    library(WithinHostExperiment); library(SummarizedExperiment) })",
        "ld <- function(f) { e <- new.env(); load(f, envir = e); e$example_whe }",
        "parts <- function(x) list(assays = as.list(assays(x)),",
        "    rowRanges = rowRanges(x), colData = colData(x),",
        "    metadata = metadata(x), qcLog = qcLog(x), qcFilters = qcFilters(x))",
        sprintf("a <- ld(%s); b <- ld(%s)", deparse(rda_before), deparse(rda)),
        "same <- identical(parts(a), parts(b))",
        "cat(if (same) 'content identical' else 'content differs',",
        "    if (identical(a, b)) '' else",
        "    '(file bytes differ: R version and class package stamps)',",
        "    '\\n')",
        "quit(status = if (same) 0L else 1L)"))
    cmp_log <- file.path(log_dir, "5_compare.log")
    st <- run_r(cmp, cmp_log, lib_env)
    same <- identical(as.integer(st), 0L)
    msg <- paste(readLines(cmp_log, warn = FALSE), collapse = " ")
    ## keep the committed bytes when only the stamps differ
    if (same) invisible(file.copy(rda_before, rda, overwrite = TRUE))
    note("5 save_example_data", if (same) "PASS" else "FAIL",
         if (same) paste(msg, "- committed file kept") else msg)
}

## ---- 6. Test suite -----------------------------------------------------
tests <- helper("run_tests.R", c(
    "suppressPackageStartupMessages({",
    "    library(testthat); library(WithinHostExperiment) })",
    "r <- as.data.frame(test_dir('tests/testthat',",
    "    package = 'WithinHostExperiment', load_package = 'installed',",
    "    reporter = 'summary', stop_on_failure = FALSE))",
    "cat('\\nTOTAL tests', nrow(r), 'expectations', sum(r$nb),",
    "    'failed', sum(r$failed), 'errors', sum(r$error),",
    "    'skipped', sum(r$skipped), '\\n')",
    "if (sum(r$failed) || any(r$error)) quit(status = 1L)"))
if (run_stage("6 tests", tests, env = lib_env, announce = FALSE)) {
    total <- grep("^TOTAL", readLines(file.path(log_dir, "6_tests.log")),
                  value = TRUE)
    note("6 tests", "PASS", trimws(sub("^TOTAL", "", total[1L])))
}

## ---- 7. README figures ------------------------------------------------
figures <- file.path("man", "figures", c("fig1_replicate_qc.png",
    "fig2_frequency_spectrum.png", "fig3_qc_overview.png",
    "fig4_diversity_landscape.png", "fig5_evolutionary_analysis.png",
    "fig6_consensus_validation.png"))
run_stage("7 generate_readme_figures",
          file.path("inst", "scripts", "generate_readme_figures.R"),
          watch = c(figures, file.path(OUT, c("readme_figure_values.csv",
                                              "figure_md5.csv",
                                              "input_md5.csv"))),
          report_only = file.path(OUT, "sessionInfo.txt"))
pdfs <- list.files(OUT, pattern = "\\.pdf$")
cat(sprintf("      vector copies in %s/: %s\n", OUT,
            if (length(pdfs)) paste(pdfs, collapse = ", ") else "none"))

## ---- 8. README text ----------------------------------------------------
run_stage("8 fill_readme", file.path("inst", "scripts", "fill_readme.R"),
          watch = "README.md")

## ---- 9. Functional checks ---------------------------------------------
run_stage("9 verify_claims", file.path("inst", "scripts", "verify_claims.R"),
          env = lib_env)

## ---- Summary -----------------------------------------------------------
cat("\n")
print(results, row.names = FALSE, right = FALSE)
failed <- results$status == "FAIL"
if (any(failed)) {
    cat("\n", sum(failed), " stage(s) failed. Logs: ", log_dir,
        "\nThe committed copies of changed files can be restored with git.\n",
        sep = "")
    quit(status = 1L)
}
cat("\nEvery stage passed: the derived files rebuild from the raw data",
    "as committed.\n")
