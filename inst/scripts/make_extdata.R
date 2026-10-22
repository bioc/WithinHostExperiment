#!/usr/bin/env Rscript
## make_extdata.R
##
## How every file in inst/extdata was made, and how to make it again.
## Run from the package root of a clone:
##   Rscript inst/scripts/make_extdata.R
##
## inst/extdata holds nineteen files in three groups.
##
## 1. Twelve synthetic fixtures, written by create_test_data.R from a seed
##    of 42 set at the top of that script: test_ref.fa (a 1000-base random
##    sequence, test_segment), test_donor.vcf, test_recipient.vcf,
##    test_independent.vcf, test_donor_rep1.vcf, test_donor_rep2.vcf,
##    test_donor_ivar.tsv, test_donor_lofreq.vcf, test_donor_freebayes.vcf,
##    test_generic.csv, test_coldata.csv and test_primers.bed. The donor,
##    recipient and independent samples carry 15, 10 and 8 iSNVs, with 5
##    shared between donor and recipient; the two replicates of the donor
##    differ at two sites and at a replicate-only site at position 500;
##    site 13 of the donor (position 812) has a depth of 150 so that the
##    depth filter removes something. Re-running the script reproduces
##    all twelve byte for byte.
##
## 2. Four excerpts of published data, written by prepare_lauring_data.R:
##    MHM1143_rep1.tsv, MHM1143_rep2.tsv, MHM1184_rep1.tsv and
##    MHM1184_rep2.tsv, the two technical replicates of the donor and the
##    recipient of household HH6 in Bendall et al. (2023), Nature
##    Communications 14:272, doi:10.1038/s41467-023-36001-5. The script
##    reads all_variants_filtered.tsv and Transmission_pairs.csv from
##    github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck, kept in
##    inst/scripts/real_data of the repository, which is not part of the
##    installed package; the upstream repository states no licence, and
##    the four files are an excerpt of 4 of its 159 calls, attributed to
##    the paper. This stage is skipped, with a message, where real_data is
##    absent.
##
## 3. Three synthetic tables with no generator: test_influenza_donor.tsv
##    (12 rows on an influenza HA segment), test_hiv_sample.tsv (20 rows on
##    an HIV-1 pol region) and test_tb_sample.tsv (3 rows on an M.
##    tuberculosis gene). They have been unchanged since the package's
##    first commit, how they were written was not recorded, and they are
##    not measurements. They are kept as they are, since regenerating them
##    would change the numbers the vignette and the tests show; this
##    script checks that they still match the checksums recorded below.
##
## data/example_whe.rda, the example object, is built from the first
## three VCFs of group 1 by save_example_data.R.
##
## inst/scripts/README.md gives the same provenance in prose.

if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] !=
        "WithinHostExperiment")
    stop("Run this script from the WithinHostExperiment package root.")

ext <- file.path("inst", "extdata")
run <- function(script) {
    cat("==", script, "\n")
    source(file.path("inst", "scripts", script), local = new.env())
}

## 1. synthetic fixtures (seed 42, inside the script)
run("create_test_data.R")

## 2. Bendall et al. excerpts, when the clone holds the inputs
if (file.exists(file.path("inst", "scripts", "real_data",
                          "all_variants_filtered.tsv"))) {
    run("prepare_lauring_data.R")
} else {
    cat("== prepare_lauring_data.R skipped: inst/scripts/real_data is not",
        "present (it is kept out of the installed package)\n")
}

## 3. the three fixtures without a generator, checked against the
##    checksums recorded here
unchanged <- c(
    test_influenza_donor.tsv = "602581c4bcb90a3e902894cb80007562",
    test_hiv_sample.tsv      = "5ef1c88c19eae2bceb04dd077f74fb1c",
    test_tb_sample.tsv       = "933d5a403149ef079b7f8fc682ff90d8")
now <- unname(tools::md5sum(file.path(ext, names(unchanged))))
for (i in seq_along(unchanged)) {
    cat(sprintf("== %-26s %s\n", names(unchanged)[i],
                if (identical(now[i], unname(unchanged[i]))) "unchanged"
                else "CHANGED: this file has no generator"))
}

## the example object, from the synthetic VCFs
run("save_example_data.R")
