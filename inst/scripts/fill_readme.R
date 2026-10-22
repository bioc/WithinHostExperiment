## fill_readme.R
##
## Writes README.md from inst/scripts/README_template.md, replacing every
## placeholder {{figure::panel::quantity::format}} with the value that
## generate_readme_figures.R recorded in
## inst/scripts/readme_values/readme_figure_values.csv. No number in the README
## is typed by hand. Run from the package root after that script:
##   Rscript inst/scripts/fill_readme.R
##
## format is a sprintf() format, "%d" (checked to be a whole number),
## "%s" (inserted as recorded), "comma" (rounded, with thousands
## separators) or "p" (a p-value as "= 0.068", or "< 0.001" when smaller
## than that, so a small p is never printed as 0.000).
## {{session::::R version::%s}} comes from sessionInfo.txt.
## The script stops if a placeholder has no recorded value, a value is not
## numeric where a number is expected, or any placeholder is left unfilled.

args <- commandArgs(trailingOnly = TRUE)
outFile <- if (length(args) >= 1L) args[[1]] else "README.md"
valuesDir <- if (length(args) >= 2L) args[[2]] else
    file.path("inst", "scripts", "readme_values")
templateFile <- file.path("inst", "scripts", "README_template.md")
valuesFile <- file.path(valuesDir, "readme_figure_values.csv")
sessionFile <- file.path(valuesDir, "sessionInfo.txt")
figureFile <- file.path(valuesDir, "figure_md5.csv")
for (f in c(templateFile, valuesFile, sessionFile, figureFile)) {
    if (!file.exists(f)) stop("missing ", f, "; run generate_readme_figures.R first")
}

values <- utils::read.csv(valuesFile, colClasses = "character", check.names = FALSE)
keys <- paste(values$figure, values$panel, values$quantity, sep = "::")
if (anyDuplicated(keys)) stop("duplicate recorded quantities: ", paste(unique(keys[duplicated(keys)]), collapse = "; "))
rVersion <- sub("^(R version [0-9.]+).*$", "\\1", readLines(sessionFile, n = 1L))
if (!grepl("^R version [0-9.]+$", rVersion)) stop("cannot read the R version from ", sessionFile)

text <- paste(readLines(templateFile, encoding = "UTF-8", warn = FALSE), collapse = "\n")
pattern <- "\\{\\{([^:{}]*)::([^:{}]*)::((?:(?!::|\\{\\{|\\}\\}).)+)::([^:{}\\s]+)\\}\\}"
hits <- gregexpr(pattern, text, perl = TRUE)
placeholders <- regmatches(text, hits)[[1]]
opened <- lengths(regmatches(text, gregexpr("\\{\\{", text)))
if (opened != length(placeholders)) {
    stop(opened - length(placeholders), " malformed placeholder(s) in the template")
}

## Every figure the README shows must be the one this run of
## generate_readme_figures.R produced, so that captions cannot be written from
## the values of one run against the figures of another.
figures <- utils::read.csv(figureFile, colClasses = "character", check.names = FALSE)
shown <- unique(regmatches(text, gregexpr('(?<=<img src=")[^"]+', text, perl = TRUE))[[1]])
for (img in shown) {
    i <- match(img, figures$file)
    if (is.na(i))
        stop("the template shows ", img, ", which is not in ", figureFile,
             "; re-run generate_readme_figures.R")
    if (!file.exists(img)) stop("missing figure ", img, "; re-run generate_readme_figures.R")
    if (!identical(unname(tools::md5sum(img)), figures$md5[i]))
        stop(img, " is not the figure recorded in ", figureFile,
             "; re-run generate_readme_figures.R so that the captions and the ",
             "figures come from one run")
}

fillOne <- function(placeholder) {
    parts <- regmatches(placeholder, regexec(pattern, placeholder, perl = TRUE))[[1]][-1]
    key <- paste(parts[1:3], collapse = "::")
    format <- parts[4]
    if (!(format %in% c("%d", "%s", "comma", "p") ||
          grepl("^%[0-9]*\\.[0-9]+[fge]$", format)))
        stop("'", key, "' asks for the format '", format,
             "', which is not one of %d, %s, comma, p or %.Nf")
    if (identical(key, "session::::R version")) return(rVersion)
    i <- match(key, keys)
    if (is.na(i)) stop("no recorded value for '", key, "'")
    value <- values$value[i]
    if (identical(format, "%s")) return(trimws(value))
    number <- suppressWarnings(as.numeric(value))
    if (is.na(number)) stop("'", key, "' is not numeric: ", value)
    if (identical(format, "comma")) return(format(round(number), big.mark = ",", scientific = FALSE))
    if (identical(format, "p")) {
        if (number < 0 || number > 1) stop("'", key, "' is not a p-value: ", value)
        return(if (number < 0.001) "< 0.001" else sprintf("= %.3f", number))
    }
    if (identical(format, "%d")) {
        if (abs(number - round(number)) > 1e-8) stop("'", key, "' is not a whole number: ", value)
        return(sprintf("%d", as.integer(round(number))))
    }
    ## A value that rounds to zero from below would otherwise print as -0.00
    sub("^-(0(\\.0*)?)$", "\\1", sprintf(format, number))
}
filled <- vapply(placeholders, fillOne, character(1), USE.NAMES = FALSE)
regmatches(text, hits) <- list(filled)
if (grepl("{{", text, fixed = TRUE) || grepl("}}", text, fixed = TRUE)) stop("unfilled placeholder remains")

con <- file(outFile, open = "wb")
writeBin(charToRaw(enc2utf8(paste0(text, "\n"))), con)
close(con)
usedKeys <- unique(vapply(placeholders, function(p) {
    paste(regmatches(p, regexec(pattern, p, perl = TRUE))[[1]][2:4], collapse = "::")
}, character(1)))
message(sprintf("Wrote %s: %d placeholders, %d distinct recorded values used, %d recorded values not quoted.",
                outFile, length(placeholders), length(setdiff(usedKeys, "session::::R version")),
                length(setdiff(keys, usedKeys))))
