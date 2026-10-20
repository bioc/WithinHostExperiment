# WithinHostExperiment 0.99.0

## Corrections made before submission

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
* Gene annotation: `annotateGenes()` labels variant sites with the
  genes they overlap, adding `GENE` to `mcols(rowRanges)`. Genes are
  supplied as a `GRanges` (a GFF3 or GTF read with
  `rtracklayer::import()`) or a `TxDb`, and the function stops, naming
  both sets, when the sites and the genes share no sequence name.
  `annotateCodonChange()` adds CDS-aware codon
  translation: given a GFF3 and reference FASTA, it determines
  the affected codon, translates both reference and alternative
  alleles using the standard genetic code, and classifies each
  SNV as Synonymous, Nonsynonymous, or Nonsense -- all in pure R
  with no external annotation tools required.
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
* Exact bottleneck estimation: `exactBottleneck()` implements
  the exact beta-binomial, Wright-Fisher, and presence-absence
  methods (Sobel Leonard et al. 2017).
* Strand bias detection: `sciStrandBias()` and
  `flagStrandBiasSCI()` for strand bias filtering that accounts
  for unequal strand coverage.
* Adaptive filtering: `adaptiveFilter()` for technology-specific
  QC presets (Illumina, ONT, PacBio HiFi).
* Diversity modelling: `fitDiversityModel()` for regression
  modelling of diversity measures, `multitestCorrection()` for
  multiple testing correction, `permutationTest()` for
  non-parametric hypothesis testing.
* Validation: `validateDiversity()` for cross-checking diversity
  estimates against literature-reported ranges.
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
