test_that("ISNVFilter constructor works with defaults", {
    f <- ISNVFilter()
    expect_s4_class(f, "ISNVFilter")
    expect_equal(f@minDepth, 100L)
    expect_equal(f@minFreq, 0.03)
    expect_equal(f@maxFreq, 0.50)
    expect_equal(f@minAltReads, 10L)
    expect_equal(f@maxStrandBias, 10.0)
    expect_identical(f@replicateConc, FALSE)
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
