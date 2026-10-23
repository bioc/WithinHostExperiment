# WithinHostExperiment 0.99.1

Changes made in response to the Bioconductor review. The package has not
been released, so functions are removed and renamed without a deprecation
cycle, and objects saved with 0.99.0 must be rebuilt: an `ISNVFilter`
saved then has a logical `replicateConc`, and a `WithinHostSFS` lacks the
`polarization` and `maxFreq` slots. The `WithinHostExperiment` validity
checks every filter in `qcFilters`, so `validObject()`, subsetting and
combining report a filter saved with 0.99.0, with a message to rebuild it
with `ISNVFilter()`, and `flagISNV()` given one names its `replicateConc`
slot.

## Interface changes

* `annotateFromGFF()` is removed, together with the package's own GFF3
  and FASTA readers. `annotateGenes(x, genes)` labels each site with the
  gene it falls in, adding a `GENE` column to `mcols(rowRanges(x))`.
  `genes` is a `GRanges`, such as a GFF3 or GTF file read with
  `rtracklayer::import()`, or a `TxDb`, from which
  `GenomicFeatures::genes()` takes the ranges; the overlap is
  `GenomicRanges::findOverlaps()`. `nameCol` selects the label column, and
  an unknown one is an error that lists those available; `collapse` keeps
  the longest overlapping gene or all of them. When the sites and the
  genes share no sequence name the function stops and prints both sets,
  so a RefSeq/GenBank naming mismatch can no longer return a column of
  `NA`.
* `annotateCodonChange()` is an S4 generic, `annotateCodonChange(x, cds,
  ref)`, taking the coding sequences as a `GRanges` of CDS records, a
  `GRangesList` or a `TxDb`, and the reference as a `DNAStringSet`; a file
  path for either is refused with a message naming the Bioconductor reader
  to use. The consequence is computed by
  `VariantAnnotation::predictCoding()`, and the output columns take its
  names: `TXNAME`, `GENEID`, `CONSEQUENCE`, `CDSLOC`, `PROTEINLOC`,
  `REFCODON`, `VARCODON`, `REFAA` and `VARAA` replace `GFF_FEATURE`,
  `REF_CODON`, `ALT_CODON`, `REF_AA`, `ALT_AA` and `AA_CLASS`. A site in
  several overlapping coding sequences is resolved by `collapse`, the
  longest transcript by default; a nucleotide that one transcript reads
  twice, at a -1 ribosomal frameshift, keeps its first reading, and a
  message says so. Three classifications change with the delegation:
  losing a stop codon is nonsynonymous rather than nonsense, an
  alternative start codon (`CTG` or `TTG`) at the first base of a coding
  sequence reads as `M`, and a site whose reference codon holds an
  ambiguous base is left unannotated, with a warning. The reference
  length is recorded in `seqlengths(rowRanges(x))`. The `GRanges` and
  `GRangesList` routes build the `TxDb` with txdbmaker; a `TxDb` built
  with `txdbmaker::makeTxDbFromGFF()` can lack transcripts the GFF3 holds
  (for NC_045512.2 it drops the ORF1ab gene, so sites in ORF1a and ORF1ab
  are left unannotated), so the man page advises passing the CDS records
  instead.
* `dndsWithinHost()` is replaced by `pNpSWithinHost(x, cds, ref)`, which
  takes the same inputs as `annotateCodonChange()` and calls it for the
  reference and variant codon of each iSNV. The quantity is a ratio of
  nonsynonymous to synonymous polymorphism within a host, not of
  substitutions between sequences, so the Jukes-Cantor correction is no
  longer applied. Nei-Gojobori sites are counted over the distinct codons
  of the genome, so codons shared by same-frame overlapping coding
  sequences, such as ORF1a and ORF1ab, are counted once, and a warning
  names the overlapping transcripts. Per gene, a shared codon counts in
  each gene that reads it, and so does an iSNV in it, each by that gene's
  own transcript. On the README's SARS-CoV-2 call set the pooled ratio
  moves from 0.2857 to 0.2867 with the correction removed. The man page
  says why the sites are counted here rather than by MSA2dist:
  `MSA2dist::codonmat2pnps()` and `dnastring2kaks()` compare aligned
  sequence pairs and cannot take per-site calls; for each of the 61 sense
  codons the sites equal those of `codonmat2pnps()`, which a test checks.
* `consensusSequences(x, ref, minFreq = 0.5)` is new. It writes each
  sample's substitutions at or above `minFreq` into the reference and
  returns one sequence per sample as a `DNAStringSet`; insertions and
  deletions are not applied, so the sequences stay aligned to the
  reference and to one another. `region` returns one coding sequence,
  joined across its parts and read in its own direction. dN/dS between
  hosts can then be estimated with `MSA2dist::dnastring2kaks()`, which is
  suggested, rather than with a model of this package's own. When `ref`
  and the sites share no sequence name it stops and names both sets, and
  a `region` on a sequence missing from `ref` is reported as such.
* `readWithinHost()` keeps its arguments but reads files with
  `VariantAnnotation::readVcf()` instead of parsing VCF text, and splits
  records with several alternative alleles with `expand()`. A per-allele
  field that the header does not declare `Number=A` on such a record stops
  the import, since `readVcf()` would give every allele the first allele's
  value; a field declared `Number=R`, whose first value is the reference
  allele's, stops it on any file, where that value was read as the
  alternative allele's count or frequency. `genome` takes a
  `DNAStringSet`, not a file path; its lengths, or those of `##contig`
  header lines, are recorded in `seqlengths()` for every sequence named,
  with a call or without. A `genome` that names none of the sequences the
  calls lie on, or is shorter than a called position, is an error, since
  every per-site statistic would divide by its length. A file holding more
  than one sample is refused; `as(vcf, "WithinHostExperiment")` converts a
  `VCF` object with all of its samples. Coercion to and from `VRanges`
  keeps every recorded sequence and its length, so the genome length
  survives a round trip, and a sequence name may contain ':', which it
  used to be split on.
* `calcDiversity()` and `compareDiversity()` lose the `genomeLength`
  argument. The length that `pi` and `watterson` divide by is read from
  `seqlengths(rowRanges(x))`, summed over every sequence recorded for
  the reference, so that a segment of a segmented genome without calls
  still counts, and a length that is not recorded is an error saying
  how to record it. `calcDiversity()` gains `n_sites_assessed`, the sites
  at which a sample has any call before QC, and `detection_limit`, the
  strictest `minFreq` of the filters applied, so that a diversity of zero
  after QC can be told apart from a sample with no calls. `example_whe`
  records its 1000-base segment length; its values are unchanged.
* `calcNeutralityTests()` loses the `genomeLength` argument it required.
  Tajima's D and Fu's Fs use the nucleotide diversity summed over the
  genome, which the frequencies give directly, so the length was divided
  in and multiplied back out and changed no output beyond rounding
  error; passing it is now an error.
* `validateDiversity()` is renamed `checkDiversityPlausibility()`. Its
  ranges are the package's own heuristics, not literature values, and the
  man page lists every one. Every sample of a `calcDiversity()` table is
  checked, where only the first row was checked before, and the result
  names the sample; a missing value is reported as not in range.
* `ISNVFilter()` loses the `minBaseQual` and `minMapQual` slots and
  arguments: the object holds no per-site quality, so they could be
  recorded but never applied, and quality filtering belongs in the
  variant caller.
* `replicateConc` holds a `ConcordanceRule`, or `NULL` for none, where
  it held a logical that was recorded but never applied.
  `ConcordanceRule(freqTolerance, requireBothDetected)` takes the two
  settings of `flagReplicateDiscordance()`, and `flagISNV()` applies it
  after the per-sample thresholds, flagging what that function flags
  and recording the same two `qcLog` rows, so the rule is stored with the
  other thresholds in `qcFilters`. A rule with no `replicate_group` of
  two or more samples to compare is an error; samples without a
  replicate are named in a warning. A logical `replicateConc` is
  refused, and so is a rule that would flag nothing: a `freqTolerance` of
  1 or more with `requireBothDetected = FALSE`, since frequencies lie in
  [0, 1]. `ConcordanceRule()` refuses a factor, string or logical
  tolerance and a `requireBothDetected` that is not `TRUE` or `FALSE`,
  rather than convert them, and `flagReplicateDiscordance()` now checks
  its arguments in the same way. `flagISNV()` checks its filter with
  `validObject()` first, so a threshold or rule changed with `@<-` into a
  value its class refuses is an error rather than a step applied.
* `consensusISNV()` keeps the `seqinfo()` of its inputs, merged, and
  stops when a sequence's length, circularity or genome differs between
  them. It takes each sample's `colData` row from the first input that
  holds the sample, and keeps a column for every input sample, `NA` where
  no site passed, with a warning naming any sample that fewer than
  `min_callers` inputs hold, since such a sample cannot reach consensus.
  Sequence lengths, host, role and pair columns used to be dropped, so
  `calcDiversity()` stopped for want of a genome length and the
  transmission functions found no pairs.
* The results of `quickBottleneck()` and `exactBottleneck()` gain
  `n_excluded`, the donor variants at or above `threshold` left out of
  the estimate, by reason (see below).

## Statistical safeguards

* `WithinHostSFS` gains a `polarization` slot: `"none"` for a folded
  spectrum, `"reference"` for an unfolded spectrum of
  alternative-allele frequencies, which treats the reference allele as
  ancestral, or `"outgroup"`. Validity ties it to `folded`, and `show()`
  prints it and labels the stored depth a median, which it is.
* `neutralityFromSFS()` stops on a folded spectrum, where it printed a
  note and returned a value, because Fay and Wu's H needs derived-allele
  frequencies; its result records the polarisation it used. It also
  stops on a spectrum without a genome length, where it returned values
  summed over the genome under the label "per site". `buildSFS()` now
  takes that length from `seqlengths(rowRanges(whe))` when it is not
  given, as `calcDiversity()` does.
* `WithinHostSFS` gains a `maxFreq` slot, the highest frequency a call
  could keep through QC. `buildSFS()` sets it to the smallest `maxFreq`
  among the filters applied with `flagISNV()` when it uses QC-passed
  calls, and to 1 otherwise. `neutralityFromSFS()` stops on an unfolded
  spectrum whose cap is below 1: the default `ISNVFilter()` removes every
  call above 0.5, which empties the classes Fay and Wu's H weights most,
  so H could not reveal an excess of high-frequency variants whatever
  the data. The spectrum must be built from the unfiltered object
  flagged with `maxFreq = 1`. `show()` prints the cap when it is below 1.
* `neutralityFromSFS()` computes pi and theta_H on the same derived-allele
  classes, `round(mid * n)` held within 1 to n - 1, with Fay and Wu's
  normalisation. Before, pi came from unrounded midpoints without the
  n/(n-1) factor, which biased H by about -theta/n under neutrality even
  on a complete spectrum. H is now centred on 0 only when every class
  from 1 to n - 1 is observed and each bin holds one class; the detection
  threshold and wide bins still shift it, to about -0.03 theta on the
  expected neutral spectrum under the defaults (20 bins, threshold 0.03,
  n = 100), much as before. pi, Tajima's D, Fu's Fs and H from this
  function change.
* `compareSFS()` stops unless both spectra share their bin breaks, their
  folding, their polarisation and their QC frequency cap, and neither is
  bias-corrected; it used to compare any two spectra with the same number
  of bins. A spectrum from before QC can therefore be compared with one
  from after it only when the filters kept `maxFreq = 1`. When either
  spectrum is empty it runs no test and returns `NA` diagnostics, where
  it returned `NaN` on a positive df.
* A site counts as segregating only when its frequency after QC lies
  strictly between 0 and 1. This applies to `n_sites`, S and theta_pi in
  `calcNeutralityTests()`, to `richness`, `watterson` and `chao1` in
  `calcDiversity()`, and to the sites `buildSFS()` bins. A frequency of 0,
  or of 1, a fixed difference that iVar reports, used to be counted
  without adding to pi, which pulled Tajima's D down.
* `tajimaD()` returns `NA`, with a warning, below three segregating
  sites, where it returned 0 at none and a number at one or two, and
  when the sample size rounds below four, where the variance of D is
  zero and it returned `Inf` or `NaN`. `fusFs()` returns `NA`, with a
  warning, when there is no segregating site or theta_pi is 0, where it
  returned 0, a number where none is defined; above n = 500 it computes
  the Ewens tail exactly, where a normal approximation made the value
  jump at the switch. `calcNeutralityTests()` reports `n_sites`, the
  number of sites each statistic rests on, documents every rule under
  which D or Fs is `NA` and warns once per cause. A sample with no calls
  gives zero diversity and zero segregating sites rather than an
  uninformative error, in `calcNeutralityTests()` and in
  `calcDiversity()`.
* The help pages of `fusFs()`, `calcNeutralityTests()` and
  `neutralityFromSFS()` say that the haplotype count is set to S + 1, the
  infinite-sites maximum. The value therefore carries no haplotype
  information and is pushed below 0 under neutrality, increasingly so as
  theta grows (-7.00 at theta = 5 with n = 100, -0.37 at theta = 1), so a
  negative Fs is not evidence of expansion.
* `pNpSWithinHost()` counts an iSNV only through a reading whose codon
  carries sites. A reading in a stop codon, including the loss of a stop
  codon, is not counted, since a stop codon carries no sites; where coding
  sequences overlap, such an iSNV is counted through a transcript that
  reads it in a sense codon, and one that every transcript reads in a stop
  codon enters neither n_S nor n_N. Each call is
  classified by translating its reference and variant codons with the
  standard genetic code and no initiator rule, so `ATG` to `CTG` or `TTG`
  at the first base of a coding sequence is nonsynonymous here, as `ATG`
  carries no synonymous site, although `annotateCodonChange()` reports it
  as synonymous. The rule applies to the totals and to each gene's
  counts. The vignette donor's ratio moves from 1.11 (3 synonymous, 11
  nonsynonymous) to 1.36 (2 and 9); the README figure values do not
  change.
* `calcSharedVariants()` gains `minDepth`, by default the strictest
  `minDepth` applied with `flagISNV()` (100 when none was), and splits the
  donor-only count into `n_donor_only_adequate_depth`,
  `n_donor_only_low_depth` and `n_donor_only_unknown_depth`, with
  `n_donor_only_recipient_failed_qc` for the sites where the recipient
  made a call that QC removed, which are not absences; only the first is
  evidence that a variant was lost in transmission. `n_donor_only` is
  still their sum.
* `exportPairFrequencies()` and `calcSharedVariants()` take a host to
  carry an allele only where its frequency is known and above 0. A site
  that a table read with `readWithinHostTable(format = "generic")` lists
  with `aaf = 0` counted as shared when no QC had been applied; it is now
  an absence, split by the recipient's depth. `recipient_freq` keeps the
  recorded 0, and sites that neither host carries are no longer listed;
  a pair with no such site gives an empty table.
* `asViralBottleneckInput()`, `quickBottleneck()` and every method of
  `exactBottleneck()` leave out, with a message, a donor variant at a site
  where QC removed a recipient call at or above `threshold`. Such a call,
  for example one above `maxFreq` where the variant rose towards fixation,
  was written as a frequency of 0 and scored as a variant lost in
  transmission, while `exact_bb` used its raw read counts. A removed call
  below `threshold` is an absence by the estimators' own definition and
  is still written as 0, and a recorded 0 is never taken for a call.
* `exactBottleneck(method = "exact_bb")` leaves out, with a message, donor
  variants whose recipient depth is unknown, which includes every site
  the recipient's VCF file does not list. They added the same constant to
  the log-likelihood at every Nb, so the estimate and the
  likelihood-ratio interval do not change, but `n_variants` counted them
  as used and a bootstrap interval resampled them; they are now counted
  in `n_excluded`, and the manual page says that `exact_bb` then sees no
  loss in transmission.
* The estimators that work from frequencies, `quickBottleneck()` and the
  `presence_absence` and `wright_fisher` methods of `exactBottleneck()`,
  use donor variants from `threshold` to `1 - threshold`, as the
  published methods do. A donor variant at frequency 1 that the recipient
  lacked gave a log-likelihood of `-Inf` at every Nb and an estimate of 1
  without a warning. Both functions now stop when the log-likelihood is
  not finite at any Nb of the grid.
* `quickBottleneck()` checks `maxNb`, where `maxNb = 0` returned an
  estimate of length 0 with two base-R warnings. It and
  `exactBottleneck()` refuse a `maxNb` that is not a single whole number
  from 1 to the integer limit; `exactBottleneck()` truncated 2.5 to 2 and
  accepted text and logical values. Both refuse a `threshold` that is not
  a single number above 0 and below 0.5: at 0 every recipient frequency
  counts as present, and from 0.5 no donor variant lies between
  `threshold` and `1 - threshold`. `quickBottleneck()` did not check it,
  and `exactBottleneck()` accepted any value from 0 to 1.
* `flagReplicateDiscordance()`, and the replicate step of `flagISNV()`,
  let a frequency difference equal to `freqTolerance` pass, as
  documented. A decimal subtraction such as 0.20 - 0.15 rounds to just
  above 0.05, so 39 of the 96 two-decimal pairs 0.05 apart were flagged;
  the comparison now allows for rounding error, and so does the frequency
  spread of `consensusISNV()`. The vignette's counts are unchanged. The
  replicate step also counts the cells it flags when `qcPass` holds
  `NA`, where the `qcLog` recorded `n_flagged` as `NA`.
* `flagTemporalInconsistency()` counts the distinct timepoints at which a
  site is detected within a host, not the samples, so technical replicates
  taken at one timepoint no longer make a variant persistent;
  `minConsecutive` runs over those timepoints, and a host sampled at a
  single timepoint is not classified. Timepoints may be numbers, dates,
  date-times, time differences, a factor in time order, or text that reads
  as numbers; other text is an error, and a blank is a missing timepoint.
  A `qcPass` cell that is `NA` counts as not passed, where an `NA` cell
  with a frequency stopped the function with an uninformative error, and
  is left `NA` when the site is flagged. Samples with no timepoint are not
  counted, with a warning.
* `flagStrandBiasSCI()` computes each sample's global strand ratio over
  the sites that sample called. The ratio was summed without removing
  missing values, so in a multi-sample object any sample with an absent
  call scored `NA` throughout and nothing in it was flagged.
  `sciStrandBias()` scores `NA` at a site with a missing count and leaves
  it out of the sums, and it refuses negative, infinite and non-numeric
  counts, where a negative count returned `-Inf`.
* `buildSFS()` checks `sampleIdx`, `nBins`, `threshold` and `fold` and
  names the one at fault; two samples are no longer pooled, and a filter
  with a `maxFreq` of 0, which leaves no call to count, is an error giving
  that value. It folds a sample with no passing
  call to an empty spectrum, where it failed and took `buildSFSList()`
  down with it. `piISNV()`, `shannonISNV()` and `simpsonISNV()` stop on
  text and on frequencies outside [0, 1], where `piISNV()` returned a
  negative diversity.
* `flagPrimerSites()` stops, naming both sets, when the primers and the
  sites share no sequence name, where it reported that no site overlapped
  a primer and flagged nothing.
* `exactBottleneck(method = "wright_fisher")` says in a message that its
  result is a heuristic score, not a likelihood; `exact_bb` and
  `presence_absence` are the estimators.

## Corrections

* `readWithinHostTable()` reads the sequence and allele columns as text.
  `read.delim()` and `read.csv()` turned a column that held only "T" into
  the logical `TRUE`, so an iVar table whose calls were all C>T was
  imported with the allele `TRUE`.
* `consensusISNV()` counts each caller once per site and sample: a caller
  that listed the same site and allele twice supplied two of the
  `min_callers` votes and was named twice in `caller_names`. It no longer
  splits a sequence name at its colons, which put every site on a contig
  such as `HLA-A*01:01` at one wrong position.
* `detectCrossContamination()` stops, listing the available columns, when
  `runCol` is not a single column of `colData`; a misspelt name used to
  compare every sample with every other without a word. Samples whose
  run is `NA` are left out with a warning giving their number.
* `plotPairScatter()` draws a donor variant whose recipient frequency is
  a recorded 0 as donor-only, where it drew it as shared, and says there
  is nothing to plot for a pair in which neither host carries an allele.

## Performance

* `consensusISNV()` builds one table of calls instead of looping over
  every cell: 3,000 sites in 20 samples from two callers took about four
  and a half minutes and take half a second, with the same result.
* `detectEmergingVariants()` orders the trajectories once and compares
  each timepoint with the next in one step, giving the same rows in the
  same order: 3,000 sites in 20 samples took about a minute and take a
  fraction of a second. Samples with no host are left out with a warning,
  where they stopped the function.

## Dependencies

* Biostrings and GenomicFeatures are imported; BiocStyle, rtracklayer,
  txdbmaker and MSA2dist are suggested. rtracklayer is not imported.
* graphics and tools, which the package already used, are declared in
  Imports; covr, which nothing used, is no longer suggested.
* `LongRead` is removed from biocViews.

## Documentation

* The vignette is rewritten as teaching material with BiocStyle: each
  code chunk is introduced by what it does and followed by how to read
  its output, terms are defined where they are introduced and collected
  in a glossary, the `qcLog` and `qcFilters` slots that
  `WithinHostExperiment` adds to `RangedSummarizedExperiment` are
  described one by one, and installation is from Bioconductor only.
* The man pages of `readWithinHost()` and `readWithinHostTable()` list
  what the import does not read. The `FILTER` column of a VCF and iVar's
  `PASS` column are not read, so records the caller rejected are imported
  as calls; per-strand counts are not read, so the strand-bias step of
  `flagISNV()` does not run on an imported object; and insertions and
  deletions are imported as width-1 rows in the caller's notation,
  without being told apart from substitutions. Each page shows how to
  drop rejected records, or keep substitutions only, before QC.
* The `WithinHostExperiment` class page states that assays must be
  ordinary R matrices, and a `qcPass` assay held as a `DelayedArray` or
  another backend is refused with a message saying so.
* The README marks installation from GitHub as not recommended and
  states the source and terms of use of the four Bendall et al. (2023)
  excerpts that ship in `inst/extdata`.
* The pkgdown reference index gives every topic a section with a
  one-line description, in place of a catch-all.
* `inst/scripts/make_extdata.R`, which ships, records how each file in
  `inst/extdata` was made and with which seed.
* `readme_figure_output/` is removed from the repository root. The values
  and checksums the README is checked against are kept in
  `inst/scripts/readme_values/`, and `inst/scripts/reproduce_all.R`
  rebuilds every derived file from the raw data and checks it against
  the committed copy.

## Tests

* `tests/testthat/test-corrections.R` holds a regression test for each
  item of 'Corrections made before submission' under 0.99.0 that still
  applies, named after the bug it guards against, and the README's
  statement that QC cannot raise pi as a test. The removed codon
  annotator is kept in `tests/testthat/helper-reference-codon.R` as a
  reference implementation, and `annotateCodonChange()` is tested
  against it.
* Each behavioural change listed above has a regression test, and most
  test names say what the behaviour used to be. The replicate rule of
  `flagISNV()` is tested at a tolerance other than the default, with
  `requireBothDetected = FALSE`, on calls that fail a threshold, at the
  tolerance boundary from both sides and on a group of three replicates,
  with hard-coded expected flags.

# WithinHostExperiment 0.99.0

## Corrections made before submission

Before submission every exported function was read against its
documentation and run on fixtures built to break it, and each
discrepancy found was fixed, documented as implemented, or both. The
items below record what changed. In 0.99.1 each item that still applies
was pinned by a regression test in `tests/testthat/test-corrections.R`,
named after the bug it guards against (for example `exact_bb matches
sites by chrom, pos and allele (was: pos only)`), so that none can
return unnoticed.

* `compareSFS()` suppressed every warning from `chisq.test()`, so a p-value
  resting on expected counts below 5 looked as trustworthy as any other. It
  now returns the smallest expected count and `chisq_approx_ok`, which is
  `FALSE` in that case, and silences only that one warning; the statistic and
  p-value are unchanged.
* `flagTemporalInconsistency()` classified each site with the last host
  processed and flagged transient sites in every sample. Sites are now
  classified within each host and flagged only in the samples of hosts where
  they are transient.
* `exactBottleneck(method = "exact_bb")` matched donor variants to rows by
  position alone, which is wrong for segmented genomes, and paired donor
  frequencies with recipient read counts by index when some sites did not
  match. Sites are now matched by chromosome, position and allele. An estimate
  that reaches `maxNb` triggers a warning.
* The documentation of `exactBottleneck()` states how `exact_bb` departs from
  Sobel Leonard et al. (2017) and that `wright_fisher` is a heuristic score.
* `flagISNV()`, `calcDiversity()` and `calcNeutralityTests()` stop on
  arguments they do not use.
* `adaptiveFilter("ont")` uses `minDepth = 60`, the depth above which Bull et
  al. (2020) found over 99% sensitivity and precision for single-nucleotide
  variants; it was 30.
* `correctSFSBias()`, `compareSFS()`, the `WithinHostSFS` bins and depth slot,
  `consensusISNV()`, `flagReplicateDiscordance()` and `validateDiversity()` are
  documented as implemented: truncation rescales the whole spectrum, bins are
  right-closed, the stored depth is a median, the consensus object holds only
  allele frequencies, and the diversity ranges are heuristic plausibility
  bounds rather than literature values.
* `ISNVFilter()` documents `minBaseQual`, `minMapQual` and `replicateConc` as
  recorded but not applied, and `flagISNV()` warns when they are set, instead
  of ignoring them silently. `WithinHostExperiment()` requires an `altFreq`
  assay, as documented.
* `trackFrequency()` builds its table in one vectorised step, which is much
  faster on genome-wide data and gives the same rows.
* `plotFrequencySpectrum()`, `plotPairScatter()` and `plotQCDashboard()` clip the
  view with coordinate limits instead of removing data outside the axis range,
  which dropped jittered points with a warning; `plotFrequencyTrajectory()`
  uses the Okabe-Ito palette, as described, and `plotSFS()` recycles it for
  more than eight samples.

## New features

* `WithinHostExperiment` S4 class extending
  `RangedSummarizedExperiment` for within-host pathogen iSNV data,
  with `setValidity` enforcing structural invariants.
* Multi-caller import: `readWithinHost()` reads VCFs in the LoFreq,
  Freebayes, generic and iVar-style field layouts with auto-detection
  (`ivar variants` writes a TSV, so an iVar-style VCF comes from a
  converter). `readWithinHostTable()` handles the iVar TSV directly,
  and generic CSV.
* Replicate-aware QC: `flagISNV()` applies configurable quality
  filters without removing data. `flagReplicateDiscordance()`
  identifies discordant technical replicates with full provenance
  logging. `flagPrimerSites()` masks amplicon primer regions.
* Diversity indices: `shannonISNV()`, `piISNV()` (with optional
  finite-sample correction per Nei 1987), `wattersonISNV()`,
  `simpsonISNV()`, `chao1ISNV()`, and `calcDiversity()` for
  per-sample within-host diversity.
* Neutrality tests: `tajimaD()`, `fusFs()`, and
  `calcNeutralityTests()` for site-frequency-spectrum-based
  neutrality testing with deep-sequencing-aware sample size
  capping.
* **Site frequency spectrum as a first-class object**:
  `WithinHostSFS` S4 class stores the within-host SFS with
  metadata (sample, genome length, depth, detection threshold).
  `buildSFS()` / `buildSFSList()` construct from a WHE;
  `correctSFSBias()` applies ascertainment bias correction for
  the detection threshold (truncation or binomial method);
  `neutralityFromSFS()` computes Tajima's D, Fu's Fs, and
  **Fay and Wu's H** (new) from a single SFS object;
  `compareSFS()` tests SFS homogeneity between samples;
  `plotSFS()` provides publication-quality SFS visualisation.
* **Temporal QC for longitudinal data**:
  `flagTemporalInconsistency()` identifies transient iSNVs that
  appear at only one timepoint (likely artefacts) versus persistent
  ones detected across multiple consecutive timepoints. Supports
  configurable `minTimepoints` and `minConsecutive` thresholds.
* **Cross-contamination detection**:
  `detectCrossContamination()` identifies index-hopping signatures
  where a high-frequency variant in one sample appears at
  suspiciously low frequency in a co-sequenced sample, consistent
  with adapter exchange during pooled sequencing. Returns a
  structured DataFrame of suspect events with source/sink
  frequencies and run metadata.
* GFF3 annotation: `annotateFromGFF()` overlaps variant sites with
  gene models from a GFF3 file, adding `GFF_FEATURE` to
  `mcols(rowRanges)` for use with `dndsWithinHost()`.
  `annotateCodonChange()` extends this with CDS-aware codon
  translation: given a GFF3 and reference FASTA, it determines
  the affected codon, translates both reference and alternative
  alleles using the standard genetic code, and classifies each
  SNV as Synonymous, Nonsynonymous, or Nonsense. (0.99.1 replaces
  both with functions that take Bioconductor objects and delegates the
  consequence call to `VariantAnnotation::predictCoding()`.)
* Selection analysis: `dndsWithinHost()` estimates within-host dN/dS
  with the Nei-Gojobori counting method, normalising by the synonymous
  and nonsynonymous sites of the CDS features in a GFF3 and reference
  FASTA, with Jukes-Cantor correction; `estimateSelectionCoefficient()`
  for time-series selection coefficient inference.
* Multi-caller consensus: `consensusISNV()` for deriving consensus
  variant calls across callers, `calcCallerConcordance()` for
  inter-caller agreement metrics, `benchmarkCallers()` for caller
  performance comparison.
* Longitudinal analysis: `trackFrequency()` for tracking variant
  allele frequencies over time, `detectEmergingVariants()` for
  identifying newly arising variants, `plotFrequencyTrajectory()`
  for publication-quality trajectory plots.
* Transmission analysis: `transmissionPairs()`,
  `exportPairFrequencies()`, `calcSharedVariants()`, and
  `quickBottleneck()` (approximate beta-binomial MLE, Sobel
  Leonard et al. 2017) for bottleneck estimation.
* Exact bottleneck estimation: `exactBottleneck()` offers a
  beta-binomial likelihood over the transmitted virions that models
  finite recipient depth, in the spirit of Sobel Leonard et al. (2017)
  with two stated departures, a presence-absence likelihood, and a
  Wright-Fisher score, which is a heuristic rather than a likelihood.
* Strand bias detection: `sciStrandBias()` and
  `flagStrandBiasSCI()` for strand bias filtering that accounts
  for unequal strand coverage.
* Adaptive filtering: `adaptiveFilter()` for technology-specific
  QC presets (Illumina, ONT, PacBio HiFi).
* Diversity modelling: `fitDiversityModel()` for regression
  modelling of diversity measures, `multitestCorrection()` for
  multiple testing correction, `permutationTest()` for
  non-parametric hypothesis testing.
* Validation: `validateDiversity()` compares diversity estimates with
  wide heuristic plausibility ranges per pathogen, not with literature
  values.
* ViralBottleneck bridge: `asViralBottleneckInput()` exports
  standardised data for the ViralBottleneck package (Zheng
  et al. 2025).
* Interoperability: `as(whe, "VRanges")`,
  `as(vr, "WithinHostExperiment")`, and
  `as(whe, "data.frame")` for seamless integration.
* Batch operations: methods for the S4Vectors `combineRows()` and
  `combineCols()` generics merge experiments together with their QC
  logs and filters.
* Publication-quality visualisation: `plotFrequencySpectrum()`,
  `plotPairScatter()`, `plotQCDashboard()`, and
  `plotFrequencyTrajectory()` with a colour-blind-safe palette
  (Wong 2011).
* SARS-CoV-2 case study in the vignette, with the same workflow run on
  synthetic examples for influenza A, HIV-1 and M. tuberculosis.
