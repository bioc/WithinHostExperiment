test_that("ISNVFilter constructor works with defaults", {
    f <- ISNVFilter()
    expect_s4_class(f, "ISNVFilter")
    expect_equal(f@minDepth, 100L)
    expect_equal(f@minFreq, 0.03)
    expect_equal(f@maxFreq, 0.50)
    expect_equal(f@minAltReads, 10L)
    expect_equal(f@maxStrandBias, 10.0)
    expect_null(f@replicateConc)
})

test_that("ConcordanceRule holds a positive tolerance and a logical", {
    r <- ConcordanceRule()
    expect_s4_class(r, "ConcordanceRule")
    expect_equal(r@freqTolerance, 0.05)
    expect_true(r@requireBothDetected)
    expect_identical(ConcordanceRule(freqTolerance = Inf)@freqTolerance, Inf)
    expect_error(ConcordanceRule(freqTolerance = 0), "freqTolerance")
    expect_error(ConcordanceRule(freqTolerance = NA), "freqTolerance")
    expect_error(ConcordanceRule(freqTolerance = c(0.05, 0.1)),
                 "freqTolerance")
    expect_error(ConcordanceRule(requireBothDetected = NA),
                 "requireBothDetected")
    ## a rule that could flag nothing is refused
    expect_error(ConcordanceRule(freqTolerance = Inf,
                                 requireBothDetected = FALSE),
                 "flags nothing")
})

test_that("ConcordanceRule refuses freqTolerance >= 1 without requireBothDetected (was: only Inf)", {
    ## frequencies lie in [0, 1], so no two differ by more than 1 and such
    ## a rule flags nothing
    expect_error(ConcordanceRule(1, FALSE), "freqTolerance >= 1")
    expect_error(ConcordanceRule(5, FALSE), "flags nothing")
    expect_error(ConcordanceRule(Inf, FALSE), "frequencies lie in \\[0, 1\\]")
    ## with requireBothDetected it is replicate intersection, and kept
    expect_s4_class(ConcordanceRule(1, TRUE), "ConcordanceRule")
    expect_s4_class(ConcordanceRule(0.99, FALSE), "ConcordanceRule")
})

test_that("ConcordanceRule() refuses a wrong type (was: coerced to a number)", {
    ## a factor became its level code and TRUE became 1; requireBothDetected
    ## took any non-zero number as TRUE
    expect_error(ConcordanceRule(factor("0.02", levels = c("0.01", "0.02"))),
                 "'freqTolerance' must be a single positive number")
    expect_error(ConcordanceRule(TRUE),
                 "'freqTolerance' must be a single positive number")
    expect_error(ConcordanceRule("0.1"),
                 "'freqTolerance' must be a single positive number")
    expect_error(ConcordanceRule(requireBothDetected = 0.5),
                 "'requireBothDetected' must be TRUE or FALSE")
    expect_error(ConcordanceRule(requireBothDetected = "T"),
                 "'requireBothDetected' must be TRUE or FALSE")
    ## an integer tolerance is a number, stored as a double
    expect_identical(ConcordanceRule(1L)@freqTolerance, 1)
})

test_that("ISNVFilter takes a ConcordanceRule or NULL as replicateConc", {
    rule <- ConcordanceRule(freqTolerance = 0.1)
    expect_identical(ISNVFilter(replicateConc = rule)@replicateConc, rule)
    ## the logical it used to record is refused, naming the class to use
    expect_error(ISNVFilter(replicateConc = TRUE),
                 "must be a ConcordanceRule or NULL")
    expect_error(ISNVFilter(replicateConc = FALSE),
                 "must be a ConcordanceRule or NULL")
    expect_error(ISNVFilter(replicateConc = 0.05),
                 "must be a ConcordanceRule or NULL")
    ## a rule changed after construction is caught by the filter's validity
    f <- ISNVFilter(replicateConc = rule)
    f@replicateConc@freqTolerance <- -1
    expect_error(validObject(f), "'replicateConc': .*freqTolerance")
})

test_that("ISNVFilter and ConcordanceRule show the replicate settings", {
    expect_output(show(ISNVFilter()), "replicateConc: none")
    expect_output(show(ISNVFilter(replicateConc = ConcordanceRule(0.1))),
                  "freqTolerance 0.1, requireBothDetected TRUE")
    expect_output(show(ConcordanceRule()), "requireBothDetected: TRUE")
    rule <- ConcordanceRule(freqTolerance = 0.25, requireBothDetected = FALSE)
    expect_output(show(rule), "freqTolerance: +0\\.25")
    expect_output(show(rule), "requireBothDetected: FALSE")
})

test_that("ISNVFilter has no base or mapping quality threshold", {
    ## those belong in the variant caller: the object holds no per-site
    ## quality for a filter to read (was: recorded but not applied)
    expect_false(any(c("minBaseQual", "minMapQual") %in%
                         slotNames("ISNVFilter")))
    expect_error(ISNVFilter(minBaseQual = 30L), "unused argument")
    expect_error(ISNVFilter(minMapQual = 30L), "unused argument")
})

test_that("ISNVFilter constructor accepts custom values", {
    f <- ISNVFilter(minDepth = 500L, minFreq = 0.05, maxFreq = 0.45)
    expect_equal(f@minDepth, 500L)
    expect_equal(f@minFreq, 0.05)
    expect_equal(f@maxFreq, 0.45)
})

test_that("ISNVFilter coerces types correctly", {
    f <- ISNVFilter(minDepth = 200, minFreq = 0.1)
    expect_true(is.integer(f@minDepth))
    expect_true(is.numeric(f@minFreq))
})

test_that("ISNVFilter rejects invalid minFreq", {
    expect_error(ISNVFilter(minFreq = -0.1), "minFreq")
    expect_error(ISNVFilter(minFreq = 1.5), "minFreq")
})

test_that("ISNVFilter rejects invalid maxFreq", {
    expect_error(ISNVFilter(maxFreq = 2.0), "maxFreq")
})

test_that("ISNVFilter rejects maxFreq < minFreq", {
    expect_error(ISNVFilter(minFreq = 0.5, maxFreq = 0.3), "maxFreq")
})

test_that("ISNVFilter rejects negative minDepth", {
    expect_error(ISNVFilter(minDepth = -10L), "minDepth")
})

test_that("ISNVFilter rejects negative minAltReads", {
    expect_error(ISNVFilter(minAltReads = -5L), "minAltReads")
})

test_that("ISNVFilter rejects non-positive maxStrandBias", {
    expect_error(ISNVFilter(maxStrandBias = 0), "maxStrandBias")
    expect_error(ISNVFilter(maxStrandBias = -1), "maxStrandBias")
})

test_that("ISNVFilter show method produces output", {
    f <- ISNVFilter()
    expect_output(show(f), "ISNVFilter object")
    expect_output(show(f), "minDepth")
    expect_output(show(f), "freq range")
})
