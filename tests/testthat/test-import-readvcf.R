# tests/testthat/test-import-readvcf.R
# Behaviour that comes with reading VCFs through VariantAnnotation::readVcf()

vcf_file <- function(lines) {
    f <- tempfile(fileext = ".vcf")
    writeLines(lines, f)
    f
}

test_that("a per-allele field not declared Number=A stops the import", {
    ## readVcf() would keep only 100 and give it to both alleles
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">",
        "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"AO\">",
        "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"RO\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\td",
        "seg\t100\t.\tA\tG,T\t100\t.\t.\tGT:AO:RO\t0/1/2:100,50:850"))
    expect_error(readWithinHost(f, colData = S4Vectors::DataFrame(
        sample_id = "d"), caller = "freebayes"),
        "AO as Number=1 rather than Number=A")
    expect_error(as(VariantAnnotation::readVcf(f), "WithinHostExperiment"),
                 "bcftools norm -m-")
})

freebayes_ao <- function(number, alt, ao) vcf_file(c("##fileformat=VCFv4.2",
    "##source=freeBayes v1.3.6",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
    "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"RO\">",
    paste0("##FORMAT=<ID=AO,Number=", number,
           ",Type=Integer,Description=\"AO\">"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\td",
    paste0("seg\t100\t.\tA\t", alt, "\t100\t.\t.\tGT:DP:RO:AO\t0/1:100:60:",
           ao)))

test_that("a FORMAT field declared Number=R stops the import (was: the reference count read as altCount)", {
    cd <- S4Vectors::DataFrame(sample_id = "d")
    ## biallelic: AO holds 60 (REF) then 40 (ALT); 60 was read as altCount
    f <- freebayes_ao("R", "G", "60,40")
    expect_error(readWithinHost(f, colData = cd, caller = "freebayes"),
                 "AO as Number=R.*Declare each such field Number=A")
    expect_error(as(VariantAnnotation::readVcf(f), "WithinHostExperiment"),
                 "the VCF declares AO as Number=R")
    ## multi-allelic: both alleles were given 30, the REF count
    f2 <- freebayes_ao("R", "G,T", "30,20,10")
    expect_error(readWithinHost(f2, colData = cd, caller = "freebayes"),
                 "AO as Number=R")
    ## Number=A, Number=1 and Number=. stay accepted on a biallelic record
    for (n in c("A", "1", "."))
        expect_identical(unname(SummarizedExperiment::assay(
            readWithinHost(freebayes_ao(n, "G", "40"), colData = cd,
                           caller = "freebayes"), "altCount")[1, 1]),
            40L, info = n)
})

test_that("an INFO field declared Number=R stops the import (was: the reference frequency read as altFreq)", {
    lofreq_af <- function(alt, af) vcf_file(c("##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=R,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        paste0("seg\t100\t.\tA\t", alt, "\t100\tPASS\tDP=1000;AF=", af)))
    cd <- S4Vectors::DataFrame(sample_id = "d")
    ## biallelic: altFreq was 0.9, the reference allele's frequency
    expect_error(readWithinHost(lofreq_af("G", "0.9,0.1"), colData = cd,
                                caller = "lofreq"), "AF as Number=R")
    ## multi-allelic: both alleles were given 0.7
    expect_error(readWithinHost(lofreq_af("G,T", "0.7,0.2,0.1"),
                                colData = cd, caller = "lofreq"),
                 "AF as Number=R")
})

test_that("Number=1 fields are fine when no record is multi-allelic", {
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    whe <- readWithinHost(f, colData = S4Vectors::DataFrame(sample_id = "d"),
                          caller = "lofreq")
    expect_equal(unname(SummarizedExperiment::assay(whe, "altFreq")[1, 1]),
                 0.1)
    expect_equal(unname(SummarizedExperiment::assay(whe, "altCount")[1, 1]),
                 100L)
})

test_that("##contig lengths are recorded as seqlengths", {
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##contig=<ID=seg,length=1000>",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    whe <- readWithinHost(f, colData = S4Vectors::DataFrame(sample_id = "d"),
                          caller = "lofreq")
    expect_identical(GenomeInfoDb::seqlengths(
        SummarizedExperiment::rowRanges(whe))[["seg"]], 1000L)
})

test_that("contig lengths that disagree across files are an error", {
    mk <- function(len) vcf_file(c("##fileformat=VCFv4.2",
        paste0("##contig=<ID=seg,length=", len, ">"),
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    expect_error(readWithinHost(c(mk(1000), mk(1200)),
        colData = S4Vectors::DataFrame(sample_id = c("a", "b")),
        caller = "lofreq"), "different lengths for seg: 1000 and 1200")
})

test_that("a genome that disagrees with a contig length is an error", {
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##contig=<ID=seg,length=1000>",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    genome <- Biostrings::DNAStringSet(c(seg = strrep("A", 900)))
    expect_error(readWithinHost(f, colData = S4Vectors::DataFrame(
        sample_id = "d"), caller = "lofreq", genome = genome),
        "seg \\(1000 vs 900\\)")
})

lofreq_two_seqs <- function() vcf_file(c("##fileformat=VCFv4.2",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "segA\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1",
    "segB\t523\t.\tC\tT\t100\tPASS\tDP=1000;AF=0.2"))

test_that("a genome naming none of the call sequences is an error (was: nothing recorded, silently)", {
    genome <- Biostrings::DNAStringSet(c(segX = strrep("A", 1000)))
    expect_error(readWithinHost(lofreq_two_seqs(),
        colData = S4Vectors::DataFrame(sample_id = "d"), caller = "lofreq",
        genome = genome), "segA, segB.*segX")
})

test_that("a genome shorter than a called position is an error (was: a wrong length, silently)", {
    genome <- Biostrings::DNAStringSet(c(segA = strrep("A", 1000),
                                         segB = strrep("A", 500)))
    expect_error(readWithinHost(lofreq_two_seqs(),
        colData = S4Vectors::DataFrame(sample_id = "d"), caller = "lofreq",
        genome = genome), "segB:523 on a sequence of 500 bases")
})

test_that("a call sequence the genome does not name is left without a length", {
    genome <- Biostrings::DNAStringSet(c(segA = strrep("A", 1000)))
    expect_warning(whe <- readWithinHost(lofreq_two_seqs(),
        colData = S4Vectors::DataFrame(sample_id = "d"), caller = "lofreq",
        genome = genome), "no sequence named segB")
    lens <- GenomeInfoDb::seqlengths(SummarizedExperiment::rowRanges(whe))
    expect_identical(lens[["segA"]], 1000L)
    expect_true(is.na(lens[["segB"]]))
    ## and calcDiversity() then says which length is missing
    expect_error(calcDiversity(whe, indices = "pi"), "not recorded for segB")
})

test_that("every recorded sequence is kept, with a call or without (was: only those with calls)", {
    ## two segments, calls on the first only
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##contig=<ID=segA,length=1000>",
        "##contig=<ID=segB,length=9000>",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "segA\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    cd <- S4Vectors::DataFrame(sample_id = "d")
    whe <- readWithinHost(f, colData = cd, caller = "lofreq")
    expect_identical(GenomeInfoDb::seqlengths(whe),
                     c(segA = 1000L, segB = 9000L))
    ## so pi divides by the whole genome, as the VRanges route does
    expect_equal(calcDiversity(whe, indices = "pi")$pi,
                 piISNV(0.1, genomeLength = 10000, meanDepth = 1000))
    ## a genome sequence without a call is recorded too
    f1 <- vcf_file(c("##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "segA\t100\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.1"))
    genome <- Biostrings::DNAStringSet(c(segA = strrep("A", 1000),
                                         segB = strrep("C", 9000)))
    whe_g <- readWithinHost(f1, colData = cd, caller = "lofreq",
                            genome = genome)
    expect_identical(GenomeInfoDb::seqlengths(whe_g),
                     c(segA = 1000L, segB = 9000L))
})

test_that("a file with more than one sample is refused, the object is not", {
    f <- vcf_file(c("##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"RO\">",
        "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"AO\">",
        "##source=freeBayes v1.3.6",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
        paste0("seg\t100\t.\tA\tG\t100\t.\tDP=2000\tGT:DP:RO:AO",
               "\t0/1:1000:900:100\t0/1:1000:700:300")))
    expect_error(readWithinHost(f, colData = S4Vectors::DataFrame(
        sample_id = "x"), caller = "freebayes"), "holds 2 samples")
    whe <- as(VariantAnnotation::readVcf(f), "WithinHostExperiment")
    expect_identical(colnames(whe), c("A", "B"))
    expect_equal(unname(SummarizedExperiment::assay(whe, "altFreq")[1, ]),
                 c(0.1, 0.3))
    expect_identical(S4Vectors::metadata(whe)$calling_info$caller,
                     "freebayes")
})

test_that("as(vcf) matches readWithinHost() on the shipped VCFs", {
    for (f in c("test_donor.vcf", "test_donor_lofreq.vcf",
                "test_donor_freebayes.vcf")) {
        path <- system.file("extdata", f, package = "WithinHostExperiment")
        from_file <- suppressMessages(readWithinHost(path,
            colData = S4Vectors::DataFrame(sample_id = "s")))
        from_obj <- as(VariantAnnotation::readVcf(path),
                       "WithinHostExperiment")
        for (a in c("altFreq", "totalDepth", "refCount", "altCount"))
            expect_identical(unname(SummarizedExperiment::assay(from_obj, a)),
                             unname(SummarizedExperiment::assay(from_file, a)),
                             info = paste(f, a))
        expect_identical(
            S4Vectors::metadata(from_obj)$calling_info$caller,
            S4Vectors::metadata(from_file)$calling_info$caller, info = f)
    }
})

test_that("a LoFreq VCF without sample columns becomes one sample", {
    path <- system.file("extdata", "test_donor_lofreq.vcf",
                        package = "WithinHostExperiment")
    whe <- as(VariantAnnotation::readVcf(path), "WithinHostExperiment")
    expect_identical(colnames(whe), "sample1")
    expect_identical(nrow(whe), 15L)
})

test_that("as(vr, 'WithinHostExperiment') keeps the sequence lengths", {
    vr <- VariantAnnotation::VRanges(seqnames = "seg",
        ranges = IRanges::IRanges(c(10, 20), width = 1),
        ref = c("A", "C"), alt = c("G", "T"),
        totalDepth = c(100L, 200L), refDepth = c(90L, 150L),
        altDepth = c(10L, 50L), sampleNames = "s1",
        seqlengths = c(seg = 500L))
    whe <- as(vr, "WithinHostExperiment")
    expect_identical(GenomeInfoDb::seqlengths(
        SummarizedExperiment::rowRanges(whe))[["seg"]], 500L)
    expect_equal(unname(SummarizedExperiment::assay(whe, "altFreq")[, 1]),
                 c(0.1, 0.25))
})

test_that("FILTER, strand counts and indels are handled as the man page says", {
    f <- vcf_file(c("##fileformat=VCFv4.2", "##source=lofreq call",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"AF\">",
        "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"DP4\">",
        "##FILTER=<ID=sb_fdr,Description=\"Strand bias\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t45\t.\tC\tT\t255\tPASS\tDP=1500;AF=0.35;DP4=488,487,262,263",
        "seg\t198\t.\tA\tG\t255\tsb_fdr\tDP=2200;AF=0.22;DP4=858,858,484,0",
        "seg\t300\t.\tA\tAGT\t255\t.\tDP=1000;AF=0.1;DP4=450,450,50,50",
        "seg\t400\t.\tCTA\tC\t255\tPASS\tDP=1000;AF=0.1;DP4=450,450,50,50"))
    whe <- flagISNV(readWithinHost(f, colData = S4Vectors::DataFrame(
        sample_id = "d"), caller = "lofreq"))
    ## the rejected record and both indels are calls that pass QC
    expect_true(all(SummarizedExperiment::assay(whe, "qcPass")))
    ## the indels are width-1 rows in the VCF's own notation
    rr <- SummarizedExperiment::rowRanges(whe)
    expect_identical(BiocGenerics::width(rr), rep(1L, 4L))
    expect_identical(rr$alt[3:4], c("AGT", "C"))
    ## no strand assays, so the strand step leaves no qcLog row
    expect_false("fwdAltCount" %in% SummarizedExperiment::assayNames(whe))
    expect_false("strand_bias_filter" %in% qcLog(whe)$step)
    expect_error(flagStrandBiasSCI(whe), "fwdAltCount")
    ## the documented FILTER recipe drops the rejected record
    vcf <- VariantAnnotation::readVcf(f)
    keep <- VariantAnnotation::filt(vcf) %in% c("PASS", ".")
    expect_identical(BiocGenerics::start(SummarizedExperiment::rowRanges(
        as(vcf[keep, ], "WithinHostExperiment"))), c(45L, 300L, 400L))
    ## and the documented subset keeps substitutions only
    sub <- whe[nchar(rr$ref) == 1L & nchar(rr$alt) == 1L, ]
    expect_identical(BiocGenerics::start(SummarizedExperiment::rowRanges(
        sub)), c(45L, 198L))
})

test_that("caller detection reads the header, given a header or a path", {
    path <- system.file("extdata", "test_donor_freebayes.vcf",
                        package = "WithinHostExperiment")
    expect_identical(WithinHostExperiment:::.detectCaller(path),
                     "freebayes")
    expect_identical(WithinHostExperiment:::.detectCaller(
        VariantAnnotation::scanVcfHeader(path)), "freebayes")
})
